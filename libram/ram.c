/* ram: read a RAM file without ROOT, a record at a time, as htslib's bam1_t.

   The record layout, the tag orders, the blocks of names, qualities and consensus, and the
   PNEXT and TLEN predictions follow ramtools' RAMNTupleRecord, RAMTags and
   RecordBlocks; a change to those needs the same change here.

   An error records why and jumps back to the public function that was called,
   which returns it. Everything allocated belongs to the RamFile, so nothing
   leaks when it does. */

#include "ram.h"
#include "rntuple.h"

#include <htscodecs/fqzcomp_qual.h>
#include <htscodecs/rANS_static4x16.h>
#include <htscodecs/tokenise_name3.h>
#include <htslib/hts_endian.h>

#include <inttypes.h>
#include <limits.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
   kFormatVersion = 9
};

/* RAMNTupleRecord's compression_flags bits. */
enum {
   kIlluminaBinning = 1 << 15,
   kDrop = 1 << 16,
   kQualInBlock = 1 << 18,
   kSeqConsensus = 1 << 20
};

/* The first byte of a quality block: rANS or fqzcomp. */
enum {
   kQualRans = 1,
   kQualFqzcomp = 2
};

/* A reference of the file that the header lacks, and RNEXT "=". */
enum {
   kAbsent = -2,
   kSame = -3
};

/* A growable byte string. */
typedef struct {
   char *s;
   size_t n, cap;
} Str;

/* A string field: its index column and its characters. */
typedef struct {
   RntCursor index, chars;
} StrField;

/* A std::vector<std::string> of METADATA's one entry. */
typedef struct {
   char **v;
   size_t n;
} Strings;

typedef struct {
   RntCursor flag, refid, pos, mapq, cigar_n, cigar, refnext, pnext_delta, tlen_delta;
   StrField seq_lit, qual;
   RntCursor seq_diff_n, seq_diff;
   RntCursor tag_int[4], tag_rg, tag_pg;
   StrField tag_str[4];
   RntCursor tag_order, tags_n, flags;
   StrField tags;
} Records;

/* The names or qualities of the current block, decoded. The block's last row
   carries them, so these cursors run ahead of the records. */
typedef struct {
   RntCursor blob_n, blob, flags; /* flags: compression_flags, for the qualities */
   uint8_t *bytes;
   char *data;
   size_t ndata, at;
   int *lengths;
   size_t nlengths, next;
} Block;

struct RamFile {
   jmp_buf trap;
   int failed;
   char err[256];
   RntFile *meta, *f;
   RntCursor mc; /* METADATA's cursors */
   StrField ms;
   sam_hdr_t *hdr;
   Strings lines, rnames, rnexts, orders, rg, pg;
   int32_t *tids, *mtids; /* rnames and rnexts as header ids */
   uint64_t *ends;
   uint64_t nblocks, entries;
   int32_t max_span;
   int sorted;
   Records r;
   Block names, quals, cons; /* cons: the consensus, whose first base cons_pos gives */
   RntCursor cons_pos;
   int64_t cons_start;
   size_t cons_len;
   /* What ram_query asked for, and where reading is. */
   int query, stopped;
   int32_t qref, qstart, qend;
   uint64_t first, block_start, block, row;
   /* The record being read. */
   Str seq, seq_lit, qual, text, tag_str[4];
   uint32_t *cigar, *diffs;
   size_t cigar_cap, diffs_cap;
   uint8_t *array, *varints;
   int64_t *values;
};

static void verror(RamFile *rf, const char *fmt, va_list ap)
{
   if (!rf->failed)
      vsnprintf(rf->err, sizeof rf->err, fmt, ap);
}

/* Records an error and returns -1, leaving the file usable. */
static int refuse(RamFile *rf, const char *fmt, ...)
{
   va_list ap;
   va_start(ap, fmt);
   verror(rf, fmt, ap);
   va_end(ap);
   return -1;
}

/* Records the first error and jumps back to the public function, which fails from then on. */
static _Noreturn void die(RamFile *rf, const char *fmt, ...)
{
   va_list ap;
   va_start(ap, fmt);
   verror(rf, fmt, ap);
   va_end(ap);
   rf->failed = 1;
   longjmp(rf->trap, 1);
}

/* Stops at the first error the RNTuple met: its cursors return zeros after one. */
static void check(RamFile *rf, const RntFile *f)
{
   const char *err = rnt_error(f);
   if (err)
      die(rf, "%s", err);
}

static void *xrealloc(RamFile *rf, void *p, size_t n)
{
   void *q = realloc(p, n ? n : 1);
   if (q == NULL)
      die(rf, "out of memory");
   return q;
}

static void reserve(RamFile *rf, Str *b, size_t extra)
{
   if (extra > SIZE_MAX / 2 - b->n)
      die(rf, "out of memory");
   if (b->n + extra > b->cap) {
      b->cap = 2 * (b->n + extra);
      b->s = xrealloc(rf, b->s, b->cap);
   }
}

static void put(RamFile *rf, Str *b, const void *s, size_t n)
{
   reserve(rf, b, n);
   memcpy(b->s + b->n, s, n);
   b->n += n;
}

static void str_field(StrField *f, RntFile *file, const char *path)
{
   rnt_cursor(&f->index, file, rnt_column(file, path, 0));
   rnt_cursor(&f->chars, file, rnt_column(file, path, 1));
}

static void str_seek(StrField *f, uint64_t row)
{
   rnt_seek(&f->index, row);
   rnt_seek_items(&f->chars, &f->index);
}

static void str_free(StrField *f)
{
   rnt_cursor_free(&f->index);
   rnt_cursor_free(&f->chars);
}

/* Reads the next string into b, replacing what it held, and ends it with '\0'. */
static void str_next(RamFile *rf, StrField *f, Str *b)
{
   const uint64_t n = rnt_len(&f->index);
   if (n > (1U << 30))
      die(rf, "string too long");
   b->n = 0;
   reserve(rf, b, (size_t)n + 1);
   rnt_bytes(&f->chars, (size_t)n, (uint8_t *)b->s);
   b->n = (size_t)n;
   b->s[b->n] = '\0';
}

static void meta_strings(RamFile *rf, const char *name, Strings *out)
{
   char path[64];
   snprintf(path, sizeof path, "%s._0", name);
   rnt_cursor_free(&rf->mc);
   str_free(&rf->ms);
   rnt_cursor(&rf->mc, rf->meta, rnt_column(rf->meta, name, 0));
   str_field(&rf->ms, rf->meta, path);
   const uint64_t n = rnt_len(&rf->mc);
   if (n > (1U << 28))
      die(rf, "bad %s", name);
   out->v = xrealloc(rf, NULL, (size_t)n * sizeof *out->v);
   for (out->n = 0; out->n < n; ++out->n) {
      str_next(rf, &rf->ms, &rf->text);
      out->v[out->n] = xrealloc(rf, NULL, rf->text.n + 1);
      memcpy(out->v[out->n], rf->text.s, rf->text.n + 1);
   }
   check(rf, rf->meta);
}

/* A number of METADATA's one entry. */
static int64_t meta_int(RamFile *rf, const char *name)
{
   rnt_cursor_free(&rf->mc);
   rnt_cursor(&rf->mc, rf->meta, rnt_column(rf->meta, name, 0));
   const int64_t v = rnt_int(&rf->mc);
   check(rf, rf->meta);
   return v;
}

static void strings_free(Strings *s)
{
   for (size_t i = 0; i < s->n; ++i)
      free(s->v[i]);
   free(s->v);
   s->v = NULL;
   s->n = 0;
}

/* The typed tags of RAMTags.h, in the order of the record's tag_ members. */
static const char *const kTypedTags[10] = {"NM", "AS", "XS", "MQ", "RG", "PG", "MD", "MC", "SA", "XA"};
enum {
   kRG = 4,
   kMD = 6
};

static void records_open(Records *r, RntFile *f)
{
   static const char *const ints[4] = {"record.tag_nm", "record.tag_as", "record.tag_xs", "record.tag_mq"};
   static const char *const strs[4] = {"record.tag_md", "record.tag_mc", "record.tag_sa", "record.tag_xa"};
   rnt_cursor(&r->flag, f, rnt_column(f, "record.flag", 0));
   rnt_cursor(&r->refid, f, rnt_column(f, "record.refid", 0));
   rnt_cursor(&r->pos, f, rnt_column(f, "record.pos", 0));
   rnt_cursor(&r->mapq, f, rnt_column(f, "record.mapq", 0));
   rnt_cursor(&r->cigar_n, f, rnt_column(f, "record.cigar", 0));
   rnt_cursor(&r->cigar, f, rnt_column(f, "record.cigar._0", 0));
   rnt_cursor(&r->refnext, f, rnt_column(f, "record.refnext", 0));
   rnt_cursor(&r->pnext_delta, f, rnt_column(f, "record.pnext_delta", 0));
   rnt_cursor(&r->tlen_delta, f, rnt_column(f, "record.tlen_delta", 0));
   str_field(&r->seq_lit, f, "record.seq_lit");
   rnt_cursor(&r->seq_diff_n, f, rnt_column(f, "record.seq_diff", 0));
   rnt_cursor(&r->seq_diff, f, rnt_column(f, "record.seq_diff._0", 0));
   str_field(&r->qual, f, "record.qual");
   for (int i = 0; i < 4; ++i) {
      rnt_cursor(&r->tag_int[i], f, rnt_column(f, ints[i], 0));
      str_field(&r->tag_str[i], f, strs[i]);
   }
   rnt_cursor(&r->tag_rg, f, rnt_column(f, "record.tag_rg", 0));
   rnt_cursor(&r->tag_pg, f, rnt_column(f, "record.tag_pg", 0));
   rnt_cursor(&r->tag_order, f, rnt_column(f, "record.tag_order", 0));
   rnt_cursor(&r->tags_n, f, rnt_column(f, "record.tags", 0));
   str_field(&r->tags, f, "record.tags._0");
   rnt_cursor(&r->flags, f, rnt_column(f, "record.compression_flags", 0));
}

static void records_seek(Records *r, uint64_t row)
{
   RntCursor *cursors[] = {&r->flag,        &r->refid,      &r->pos,        &r->mapq,       &r->refnext,
                           &r->pnext_delta, &r->tlen_delta, &r->tag_int[0], &r->tag_int[1], &r->tag_int[2],
                           &r->tag_int[3],  &r->tag_rg,     &r->tag_pg,     &r->tag_order,  &r->flags};
   for (size_t i = 0; i < sizeof cursors / sizeof *cursors; ++i)
      rnt_seek(cursors[i], row);
   rnt_seek(&r->cigar_n, row);
   rnt_seek_items(&r->cigar, &r->cigar_n);
   rnt_seek(&r->seq_diff_n, row);
   rnt_seek_items(&r->seq_diff, &r->seq_diff_n);
   StrField *strs[] = {&r->seq_lit, &r->qual, &r->tag_str[0], &r->tag_str[1], &r->tag_str[2], &r->tag_str[3]};
   for (size_t i = 0; i < sizeof strs / sizeof *strs; ++i)
      str_seek(strs[i], row);
   rnt_seek(&r->tags_n, row);
   rnt_seek_items(&r->tags.index, &r->tags_n);
   rnt_seek_items(&r->tags.chars, &r->tags.index);
}

static void records_free(Records *r)
{
   RntCursor *cursors[] = {&r->flag,       &r->refid,      &r->pos,         &r->mapq,       &r->cigar_n,
                           &r->cigar,      &r->refnext,    &r->pnext_delta, &r->tlen_delta, &r->tag_int[0],
                           &r->tag_int[1], &r->tag_int[2], &r->tag_int[3],  &r->tag_rg,     &r->tag_pg,
                           &r->tag_order,  &r->tags_n,     &r->flags,       &r->seq_diff_n, &r->seq_diff};
   for (size_t i = 0; i < sizeof cursors / sizeof *cursors; ++i)
      rnt_cursor_free(cursors[i]);
   StrField *strs[] = {&r->seq_lit, &r->qual, &r->tag_str[0], &r->tag_str[1], &r->tag_str[2], &r->tag_str[3], &r->tags};
   for (size_t i = 0; i < sizeof strs / sizeof *strs; ++i)
      str_free(strs[i]);
}

/* Reads the blob of a block of rows, skipping the empty ones before it. */
static size_t read_blob(RamFile *rf, Block *b, uint64_t rows)
{
   size_t n = 0;
   for (uint64_t i = 0; i < rows; ++i) {
      const uint64_t len = rnt_len(&b->blob_n);
      if (len > (1U << 30))
         die(rf, "block too large");
      if (i + 1 < rows) {
         rnt_bytes(&b->blob, (size_t)len, NULL);
         continue;
      }
      n = (size_t)len;
      b->bytes = xrealloc(rf, b->bytes, n);
      rnt_bytes(&b->blob, n, b->bytes);
   }
   check(rf, rf->f);
   return n;
}

/* NameBlockReader: the block's names, each ended by '\0'. */
static void load_names(RamFile *rf, uint64_t rows)
{
   Block *b = &rf->names;
   if (b->data && b->at != b->ndata)
      die(rf, "a block has more names than records");
   const size_t n = read_blob(rf, b, rows);
   uint32_t size = 0;
   free(b->data);
   b->data = (char *)tok3_decode_names(b->bytes, (uint32_t)n, &size);
   if (b->data == NULL)
      die(rf, "the name tokeniser could not decode a block");
   uint64_t names = 0;
   for (uint32_t i = 0; i < size; ++i)
      names += b->data[i] == '\0';
   if (names != rows || (size > 0 && b->data[size - 1] != '\0'))
      die(rf, "a block has %" PRIu64 " names for %" PRIu64 " records", names, rows);
   b->ndata = size;
   b->at = 0;
}

static uint8_t *unrans(RamFile *rf, uint8_t *in, size_t n, unsigned *out)
{
   uint8_t *p = rans_uncompress_4x16(in, (unsigned)n, out);
   if (p == NULL)
      die(rf, "rANS could not decode a quality block");
   return p;
}

/* QualityBlockReader: the lossless qualities, without the SAM offset, of the block's records that have
   kQualInBlock. */
static void load_quals(RamFile *rf, uint64_t rows)
{
   Block *b = &rf->quals;
   if (b->next != b->nlengths)
      die(rf, "a block has more qualities than records");
   size_t records = 0;
   for (uint64_t i = 0; i < rows; ++i)
      records += (rnt_int(&b->flags) & kQualInBlock) != 0;
   const size_t n = read_blob(rf, b, rows);
   free(b->data);
   b->data = NULL;
   b->ndata = b->at = b->next = 0;
   b->nlengths = records;
   b->lengths = xrealloc(rf, b->lengths, records * sizeof *b->lengths);
   if (records == 0)
      return;
   if (n < 2)
      die(rf, "a quality block is damaged");
   uint8_t *in = b->bytes;
   size_t total = 0;
   if (in[0] == kQualRans) {
      /* The size of the lengths part, the lengths as varints, then the qualities: each part rANS. */
      if (n < 5)
         die(rf, "a quality block is damaged");
      const size_t nlen = (size_t)in[1] | (size_t)in[2] << 8 | (size_t)in[3] << 16 | (size_t)in[4] << 24;
      if (nlen >= n - 5)
         die(rf, "a quality block is damaged");
      unsigned nv = 0;
      free(rf->varints);
      rf->varints = unrans(rf, in + 5, nlen, &nv);
      const uint8_t *v = rf->varints;
      size_t k = 0;
      uint32_t x = 0;
      int shift = 0;
      for (unsigned i = 0; i < nv; ++i) {
         if (shift > 28)
            die(rf, "a quality block is damaged");
         x |= (uint32_t)(v[i] & 0x7f) << shift;
         shift += 7;
         if (v[i] < 0x80) {
            if (k == records || x > (uint32_t)INT32_MAX)
               die(rf, "a quality block is damaged");
            b->lengths[k++] = (int)x;
            total += x;
            x = 0;
            shift = 0;
         }
      }
      unsigned nq = 0;
      b->data = (char *)unrans(rf, in + 5 + nlen, n - 5 - nlen, &nq);
      b->ndata = nq;
      if (shift != 0 || k != records || nq != total)
         die(rf, "a quality block is damaged");
   } else if (in[0] == kQualFqzcomp) {
      size_t size = 0;
      b->data = fqz_decompress((char *)in + 1, n - 1, &size, b->lengths, (int)records);
      if (b->data == NULL)
         die(rf, "fqzcomp could not decode a quality block");
      for (size_t i = 0; i < records; ++i) {
         if (b->lengths[i] < 0)
            die(rf, "a quality block is damaged");
         total += (size_t)b->lengths[i];
      }
      if (total != size)
         die(rf, "a quality block is damaged");
      b->ndata = size;
   } else {
      die(rf, "a quality block has an unknown codec");
   }
}

/* SeqBlockReader: the block's consensus and the position of its first base. */
static void load_consensus(RamFile *rf, uint64_t rows)
{
   rf->cons_len = read_blob(rf, &rf->cons, rows);
   for (uint64_t i = 0; i < rows; ++i)
      rf->cons_start = rnt_int(&rf->cons_pos);
   check(rf, rf->f);
}

/* Rebuilds into rf->seq the SEQ that a record with kSeqConsensus codes against its block's consensus, as
   RecordBlocks' DecodeSeq() does. */
static void decode_seq(RamFile *rf, const char *name, int32_t pos, size_t ncigar, size_t ndiffs)
{
   static const char kAcgt[] = "ACGT";
   static const char kBases[] = "=ACMGRSVTWYHKDBN";
   /* Each difference is its step from the previous one, shifted left 5 bits, ORed with its code. */
   enum {
      kSubOther = 3,
      kDiffCodeBits = 5
   };
   Str *s = &rf->seq;
   const Str *lit = &rf->seq_lit;
   const uint8_t *cons = rf->cons.bytes;
   int64_t at = (int64_t)pos - rf->cons_start;
   size_t used = 0;
   s->n = 0;
   for (size_t i = 0; i < ncigar; ++i) {
      const int type = bam_cigar_type(bam_cigar_op(rf->cigar[i]));
      const size_t n = bam_cigar_oplen(rf->cigar[i]);
      if (type == 3) {
         if (at < 0 || (uint64_t)at > rf->cons_len || n > rf->cons_len - (size_t)at)
            die(rf, "the SEQ of record %s is damaged", name);
         put(rf, s, cons + at, n);
      } else if (type & 1) {
         if (n > lit->n - used)
            die(rf, "the SEQ of record %s is damaged", name);
         put(rf, s, lit->s + used, n);
         used += n;
      }
      if (type & 2)
         at += (int64_t)n;
   }
   if (used != lit->n)
      die(rf, "the SEQ of record %s is damaged", name);
   size_t p = 0;
   for (size_t k = 0; k < ndiffs; ++k) {
      p += rf->diffs[k] >> kDiffCodeBits;
      if (p >= s->n)
         die(rf, "the SEQ of record %s is damaged", name);
      const unsigned code = rf->diffs[k] & ((1U << kDiffCodeBits) - 1);
      if (code < kSubOther) {
         const char *c = s->s[p] ? strchr(kAcgt, s->s[p]) : NULL;
         if (c == NULL)
            die(rf, "the SEQ of record %s is damaged", name);
         s->s[p] = kAcgt[((size_t)(c - kAcgt) + code + 1) % 4];
      } else if (code - kSubOther < 16) {
         s->s[p] = kBases[code - kSubOther];
      } else {
         die(rf, "the SEQ of record %s is damaged", name);
      }
   }
   reserve(rf, s, 1);
   s->s[s->n] = '\0';
}

static void block_free(Block *b)
{
   rnt_cursor_free(&b->blob_n);
   rnt_cursor_free(&b->blob);
   rnt_cursor_free(&b->flags);
   free(b->bytes);
   free(b->data);
   free(b->lengths);
}

static int32_t find_ref(const Strings *rnames, const char *name, size_t len)
{
   for (size_t i = 0; i < rnames->n && i <= INT32_MAX; ++i) {
      if (strlen(rnames->v[i]) == len && memcmp(rnames->v[i], name, len) == 0)
         return (int32_t)i;
   }
   return -1;
}

/* A position in a region: digits only, at most INT32_MAX; -1 otherwise. */
static int64_t position(const char *s, size_t n)
{
   int64_t v = 0;
   for (size_t i = 0; i < n; ++i) {
      if (s[i] < '0' || s[i] > '9' || v > INT32_MAX)
         return -1;
      v = v * 10 + (s[i] - '0');
   }
   return n > 0 && v <= INT32_MAX ? v : -1;
}

/* ramntuplescan()'s reading of a region: the whole of it as a reference name, else rname[:start[-end]],
   1-based and inclusive. Start and end come back 0-based; -1 for a region the file cannot have. */
static int32_t parse_region(RamFile *rf, const char *region, int32_t *start, int32_t *end)
{
   *start = 0;
   *end = INT32_MAX - 1;
   int32_t refid = find_ref(&rf->rnames, region, strlen(region));
   if (refid >= 0)
      return refid;
   const char *colon = strrchr(region, ':');
   const size_t name_len = colon ? (size_t)(colon - region) : strlen(region);
   if (colon) {
      const char *dash = strchr(colon + 1, '-');
      const int64_t a = position(colon + 1, dash ? (size_t)(dash - colon - 1) : strlen(colon + 1));
      const int64_t b = dash ? position(dash + 1, strlen(dash + 1)) : a;
      if (a < 0 || b < 0)
         return refuse(rf, "invalid region %s: use rname[:start[-end]]", region);
      *start = a > 0 ? (int32_t)(a - 1) : 0;
      *end = (int32_t)(b - 1);
   }
   refid = find_ref(&rf->rnames, region, name_len);
   if (refid < 0)
      return refuse(rf, "reference %.*s not found", (int)name_len, region);
   return refid;
}

/* The first row at or after (refid, pos) in a sorted file, where unplaced records sort last. */
static uint64_t first_row_at_or_after(RntFile *f, uint64_t total, int32_t refid, int32_t pos)
{
   RntCursor rc;
   RntCursor pc;
   rnt_cursor(&rc, f, rnt_column(f, "record.refid", 0));
   rnt_cursor(&pc, f, rnt_column(f, "record.pos", 0));
   uint64_t lo = 0;
   uint64_t hi = total;
   while (lo < hi) {
      const uint64_t mid = lo + (hi - lo) / 2;
      rnt_seek(&rc, mid);
      const int64_t r = rnt_int(&rc);
      int before = r >= 0 && r < refid;
      if (r == refid) {
         rnt_seek(&pc, mid);
         before = rnt_int(&pc) < pos;
      }
      if (before)
         lo = mid + 1;
      else
         hi = mid;
   }
   rnt_cursor_free(&rc);
   rnt_cursor_free(&pc);
   return lo;
}

/* RAMNTupleRecord::GetRefSpan(): the reference bases the CIGAR covers. */
static int64_t ref_span(const uint32_t *cigar, size_t n)
{
   uint32_t span = 0;
   for (size_t i = 0; i < n; ++i) {
      const uint32_t op = cigar[i] & 0xf;
      if (op == 0 || op == 2 || op == 3 || op == 7 || op == 8)
         span += cigar[i] >> 4;
   }
   return span;
}

/* Modulo 2^32, as RAMNTupleRecord's Wrap(). */
static int32_t wrap(int64_t v)
{
   return (int32_t)(uint32_t)(uint64_t)v;
}

/* RAMNTupleRecord::PredictTLEN(). */
static int32_t predict_tlen(uint16_t flag, int32_t pos, int64_t span, int32_t pnext0, int same_ref)
{
   if (!same_ref || (flag & 0x4) || (flag & 0x8) || span == 0)
      return 0;
   const int64_t own5 = (flag & 0x10) ? pos + span - 1 : pos;
   const int64_t mate5 = (flag & 0x20) ? pnext0 + span - 1 : pnext0;
   const int64_t distance = mate5 - own5;
   if (distance > 0)
      return wrap(distance + 1);
   if (distance < 0)
      return wrap(distance - 1);
   return 0;
}

/* The header's id for a reference of the file: -1 for "*" or an id the file does not have. */
static int32_t tid_of(RamFile *rf, const int32_t *tids, const Strings *names, int64_t id)
{
   if (id < 0 || (uint64_t)id >= names->n)
      return -1;
   if (tids[id] == kAbsent)
      die(rf, "reference %s is not in the header", names->v[id]);
   return tids[id];
}

static int32_t *header_ids(RamFile *rf, const Strings *names)
{
   int32_t *ids = xrealloc(rf, NULL, names->n * sizeof *ids);
   for (size_t i = 0; i < names->n; ++i) {
      const int tid = sam_hdr_name2tid(rf->hdr, names->v[i]);
      if (strcmp(names->v[i], "*") == 0)
         ids[i] = -1;
      else if (strcmp(names->v[i], "=") == 0)
         ids[i] = kSame;
      else
         ids[i] = tid >= 0 ? tid : kAbsent;
   }
   return ids;
}

static void aux(RamFile *rf, bam1_t *b, const char *tag, char type, size_t len, const void *data)
{
   if (len > INT_MAX || bam_aux_append(b, tag, type, (int)len, data) < 0)
      die(rf, "cannot add tag %.2s to a record", tag);
}

/* An integer tag in the smallest type that holds it, as sam_parse1 stores one. */
static void int_tag(RamFile *rf, bam1_t *b, const char *tag, int64_t v)
{
   if (v < INT32_MIN || v > UINT32_MAX)
      die(rf, "tag %.2s is out of range", tag);
   static const char kTypes[2][3] = {{'C', 'S', 'I'}, {'c', 's', 'i'}};
   const size_t len = v >= INT8_MIN && v <= UINT8_MAX ? 1 : v >= INT16_MIN && v <= UINT16_MAX ? 2 : 4;
   uint8_t le[4];
   u32_to_le((uint32_t)v, le);
   aux(rf, b, tag, kTypes[v < 0][len / 2], len, le);
}

/* A number of a text tag, from s up to end, that must fit [lo, hi]. */
static int64_t tag_number(RamFile *rf, const char *s, const char *end, int64_t lo, int64_t hi)
{
   char *stop = NULL;
   const long long v = s < end ? strtoll(s, &stop, 10) : 0;
   if (stop != end || v < lo || v > hi)
      die(rf, "bad number in a tag: %.*s", (int)(end - s), s);
   return v;
}

static float tag_float(RamFile *rf, const char *s, const char *end)
{
   char *stop = NULL;
   const float v = s < end ? strtof(s, &stop) : 0;
   if (stop != end)
      die(rf, "bad number in a tag: %.*s", (int)(end - s), s);
   return v;
}

/* B:t,v,...: the array as BAM stores it, its type, its count and its values. Integers out of the type's range
   widen it to the smallest type that holds them all, as sam_parse1 does. */
static void array_tag(RamFile *rf, bam1_t *b, const char *tag, const char *s, const char *end)
{
   static const char kTypes[] = "cCsSiIf";
   static const int64_t kLo[] = {INT8_MIN, 0, INT16_MIN, 0, INT32_MIN, 0};
   static const int64_t kHi[] = {INT8_MAX, UINT8_MAX, INT16_MAX, UINT16_MAX, INT32_MAX, UINT32_MAX};
   const char *t = s < end && *s != '\0' ? strchr(kTypes, *s) : NULL;
   if (t == NULL)
      die(rf, "bad array type in tag %.2s", tag);
   size_t type = (size_t)(t - kTypes);
   size_t count = 0;
   for (const char *p = s + 1; p < end; ++p)
      count += *p == ',';
   if (count > UINT32_MAX || count > (SIZE_MAX - 5) / sizeof *rf->values)
      die(rf, "array too long in tag %.2s", tag);
   rf->array = xrealloc(rf, rf->array, 5 + count * 4);
   rf->values = xrealloc(rf, rf->values, count * sizeof *rf->values);
   int64_t lo = 0;
   int64_t hi = 0;
   const char *p = s + 1;
   for (size_t i = 0; i < count; ++i) {
      if (*p != ',')
         die(rf, "bad array in tag %.2s", tag);
      const char *next = memchr(p + 1, ',', (size_t)(end - p - 1));
      if (next == NULL)
         next = end;
      if (type == 6) {
         float_to_le(tag_float(rf, p + 1, next), rf->array + 5 + 4 * i);
      } else {
         const int64_t v = rf->values[i] = tag_number(rf, p + 1, next, INT32_MIN, UINT32_MAX);
         lo = i == 0 || v < lo ? v : lo;
         hi = i == 0 || v > hi ? v : hi;
      }
      p = next;
   }
   if (p != end)
      die(rf, "bad array in tag %.2s", tag);
   if (type < 6 && (lo < kLo[type] || hi > kHi[type])) {
      for (type = lo < 0 ? 0 : 1; type < 6 && (lo < kLo[type] || hi > kHi[type]); type += 2) {
      }
      if (type == 6)
         die(rf, "array out of range in tag %.2s", tag);
   }
   const size_t size = type < 2 ? 1 : type < 4 ? 2 : 4;
   for (size_t i = 0; type < 6 && i < count; ++i) {
      uint8_t le[4];
      u32_to_le((uint32_t)rf->values[i], le);
      memcpy(rf->array + 5 + i * size, le, size);
   }
   rf->array[0] = (uint8_t)kTypes[type];
   u32_to_le((uint32_t)count, rf->array + 1);
   aux(rf, b, tag, 'B', 5 + count * size, rf->array);
}

/* A tag kept as SAM text, TG:t:value, ended by '\0', stored as sam_parse1 stores it. */
static void text_tag(RamFile *rf, bam1_t *b, const char *s, size_t n)
{
   if (n < 5 || s[2] != ':' || s[4] != ':')
      die(rf, "bad tag %s", s);
   const char *v = s + 5;
   const char *end = s + n;
   uint8_t le[4];
   switch (s[3]) {
   case 'A':
      if (n != 6)
         die(rf, "bad tag %s", s);
      aux(rf, b, s, 'A', 1, v);
      break;
   case 'i': int_tag(rf, b, s, tag_number(rf, v, end, INT32_MIN, UINT32_MAX)); break;
   case 'f':
      float_to_le(tag_float(rf, v, end), le);
      aux(rf, b, s, 'f', 4, le);
      break;
   case 'Z':
   case 'H': aux(rf, b, s, s[3], n - 5 + 1, v); break;
   case 'B': array_tag(rf, b, s, v, end); break;
   default: die(rf, "bad tag %s", s);
   }
}

/* Reads from the start of the block that holds the first row ram_query asked for. The cursors move there
   when the first record is read, so a query that replaces this one costs nothing. */
static void start(RamFile *rf, int query, int32_t qref, int32_t qstart, int32_t qend)
{
   rf->query = query;
   rf->qref = qref;
   rf->qstart = qstart;
   rf->qend = qend;
   rf->stopped = 0;
   /* In a sorted file the scan starts a read's longest span before the region; otherwise it reads every
      record. */
   const int32_t from = rf->max_span > 0 && qstart > rf->max_span ? qstart - rf->max_span : 0;
   rf->first = query && rf->sorted ? first_row_at_or_after(rf->f, rf->entries, qref, from) : 0;
   uint64_t block = 0;
   while (block < rf->nblocks && rf->ends[block] < rf->first)
      block++;
   rf->block = block;
   rf->block_start = rf->row = block == 0 ? 0 : rf->ends[block - 1] + 1;
   free(rf->names.data);
   rf->names.data = NULL;
   rf->names.ndata = rf->names.at = 0;
   rf->quals.nlengths = rf->quals.next = 0;
   check(rf, rf->f);
}

static void open_file(RamFile *rf, const char *path)
{
   char err[256];
   rf->meta = rnt_open(path, "METADATA", err, sizeof err);
   if (rf->meta == NULL)
      die(rf, "%s: %s", path, err);
   const int64_t version = rnt_entries(rf->meta) ? meta_int(rf, "format_version") : -1;
   if (version != kFormatVersion)
      die(rf, "%s has RAM format version %" PRId64 "; this reader reads version %d", path, version, kFormatVersion);
   meta_strings(rf, "header_lines", &rf->lines);
   meta_strings(rf, "rname_refs", &rf->rnames);
   meta_strings(rf, "rnext_refs", &rf->rnexts);
   meta_strings(rf, "tag_orders", &rf->orders);
   meta_strings(rf, "tag_rg", &rf->rg);
   meta_strings(rf, "tag_pg", &rf->pg);
   rf->max_span = wrap(meta_int(rf, "max_ref_span"));
   rf->sorted = meta_int(rf, "coordinate_sorted") != 0;
   rnt_cursor_free(&rf->mc);
   str_free(&rf->ms);
   rnt_cursor(&rf->mc, rf->meta, rnt_column(rf->meta, "block_ends", 0));
   rf->nblocks = rnt_len(&rf->mc);
   if (rf->nblocks > (1U << 28))
      die(rf, "bad block_ends");
   rf->ends = xrealloc(rf, NULL, (size_t)rf->nblocks * sizeof *rf->ends);
   rnt_cursor(&rf->ms.index, rf->meta, rnt_column(rf->meta, "block_ends._0", 0));
   for (uint64_t i = 0; i < rf->nblocks; ++i)
      rf->ends[i] = (uint64_t)rnt_int(&rf->ms.index);
   check(rf, rf->meta);
   rnt_cursor_free(&rf->mc);
   rnt_cursor_free(&rf->ms.index);
   rnt_close(rf->meta);
   rf->meta = NULL;

   /* The header as the input had it, and the file's references as its ids. */
   rf->text.n = 0;
   for (size_t i = 0; i < rf->lines.n; ++i) {
      put(rf, &rf->text, rf->lines.v[i], strlen(rf->lines.v[i]));
      put(rf, &rf->text, "\n", 1);
   }
   put(rf, &rf->text, "", 1);
   rf->hdr = sam_hdr_parse(rf->text.n - 1, rf->text.s);
   if (rf->hdr == NULL)
      die(rf, "%s: the header cannot be parsed", path);
   strings_free(&rf->lines);
   rf->tids = header_ids(rf, &rf->rnames);
   rf->mtids = header_ids(rf, &rf->rnexts);

   rf->f = rnt_open(path, "RAM", err, sizeof err);
   if (rf->f == NULL)
      die(rf, "%s: %s", path, err);
   rf->entries = rnt_entries(rf->f);
   for (uint64_t i = 0; i < rf->nblocks; ++i) {
      if ((i > 0 && rf->ends[i] <= rf->ends[i - 1]) || rf->ends[i] >= rf->entries)
         die(rf, "the block ends in the file are damaged");
   }
   if (rf->entries > 0 && (rf->nblocks == 0 || rf->ends[rf->nblocks - 1] != rf->entries - 1))
      die(rf, "records after the last block");
   records_open(&rf->r, rf->f);
   rnt_cursor(&rf->names.blob_n, rf->f, rnt_column(rf->f, "nameblock", 0));
   rnt_cursor(&rf->names.blob, rf->f, rnt_column(rf->f, "nameblock._0", 0));
   rnt_cursor(&rf->quals.blob_n, rf->f, rnt_column(rf->f, "qualblock", 0));
   rnt_cursor(&rf->quals.blob, rf->f, rnt_column(rf->f, "qualblock._0", 0));
   rnt_cursor(&rf->quals.flags, rf->f, rnt_column(rf->f, "record.compression_flags", 0));
   rnt_cursor(&rf->cons.blob_n, rf->f, rnt_column(rf->f, "consensus", 0));
   rnt_cursor(&rf->cons.blob, rf->f, rnt_column(rf->f, "consensus", 1));
   rnt_cursor(&rf->cons_pos, rf->f, rnt_column(rf->f, "consensus_pos", 0));
   check(rf, rf->f);
   start(rf, 0, -1, 0, 0);
}

static int next_record(RamFile *rf, bam1_t *b)
{
   Records *r = &rf->r;
   Block *names = &rf->names;
   Block *quals = &rf->quals;
   while (rf->row < rf->entries && !rf->stopped) {
      const uint64_t row = rf->row++;
      if (row == rf->block_start) {
         records_seek(r, row);
         rnt_seek(&names->blob_n, row);
         rnt_seek_items(&names->blob, &names->blob_n);
         rnt_seek(&quals->blob_n, row);
         rnt_seek_items(&quals->blob, &quals->blob_n);
         rnt_seek(&quals->flags, row);
         rnt_seek(&rf->cons.blob_n, row);
         rnt_seek_items(&rf->cons.blob, &rf->cons.blob_n);
         rnt_seek(&rf->cons_pos, row);
      }
      if (row == rf->block_start || row == rf->ends[rf->block - 1] + 1) {
         load_names(rf, rf->ends[rf->block] - row + 1);
         load_quals(rf, rf->ends[rf->block] - row + 1);
         load_consensus(rf, rf->ends[rf->block] - row + 1);
         rf->block++;
      }
      const uint16_t flag = (uint16_t)rnt_int(&r->flag);
      const int64_t refid = rnt_int(&r->refid);
      const int32_t pos = (int32_t)rnt_int(&r->pos);
      const int64_t mapq = rnt_int(&r->mapq);
      const uint64_t ncigar = rnt_len(&r->cigar_n);
      if (ncigar > (1U << 28))
         die(rf, "CIGAR too long");
      if (ncigar > rf->cigar_cap) {
         rf->cigar_cap = (size_t)ncigar;
         rf->cigar = xrealloc(rf, rf->cigar, rf->cigar_cap * sizeof *rf->cigar);
      }
      for (size_t i = 0; i < ncigar; ++i)
         rf->cigar[i] = (uint32_t)rnt_int(&r->cigar);
      const int64_t refnext = rnt_int(&r->refnext);
      const int32_t pnext_delta = (int32_t)rnt_int(&r->pnext_delta);
      const int32_t tlen_delta = (int32_t)rnt_int(&r->tlen_delta);
      str_next(rf, &r->seq_lit, &rf->seq_lit);
      const uint64_t ndiffs = rnt_len(&r->seq_diff_n);
      if (ndiffs > (1U << 30))
         die(rf, "too many SEQ differences");
      if (ndiffs > rf->diffs_cap) {
         rf->diffs_cap = (size_t)ndiffs;
         rf->diffs = xrealloc(rf, rf->diffs, rf->diffs_cap * sizeof *rf->diffs);
      }
      for (size_t i = 0; i < ndiffs; ++i)
         rf->diffs[i] = (uint32_t)rnt_int(&r->seq_diff);
      str_next(rf, &r->qual, &rf->qual);
      int32_t tag_int[4];
      for (int i = 0; i < 4; ++i)
         tag_int[i] = (int32_t)rnt_int(&r->tag_int[i]);
      const int64_t tag_rg = rnt_int(&r->tag_rg);
      const int64_t tag_pg = rnt_int(&r->tag_pg);
      for (int i = 0; i < 4; ++i)
         str_next(rf, &r->tag_str[i], &rf->tag_str[i]);
      const int64_t tag_order = rnt_int(&r->tag_order);
      const uint64_t ntags = rnt_len(&r->tags_n);
      const uint32_t flags = (uint32_t)rnt_int(&r->flags);
      check(rf, rf->f);

      /* Whether ramntuplescan() keeps the row. The name, the qualities and the text tags are read either
         way: the blocks and the cursors move one row at a time. */
      const int64_t span = ref_span(rf->cigar, (size_t)ncigar);
      int keep = 1;
      if (rf->query) {
         /* A record with no reference span still covers the base it sits on, as htslib's bam_endpos says. */
         keep = refid == rf->qref && pos <= rf->qend &&
                (pos >= rf->qstart || pos + (span > 0 ? span : 1) - 1 >= rf->qstart);
         if (rf->sorted && row >= rf->first && (refid < 0 || refid > rf->qref || (refid == rf->qref && pos > rf->qend)))
            rf->stopped = 1;
         keep = keep && row >= rf->first && !rf->stopped;
      }
      const char *name = names->data + names->at;
      const size_t name_len = strlen(name);
      names->at += name_len + 1;
      const char *block_qual = NULL;
      size_t block_qual_len = 0;
      if (flags & kQualInBlock) {
         if (quals->next == quals->nlengths)
            die(rf, "more qualities in records than in their block");
         block_qual_len = (size_t)quals->lengths[quals->next++];
         block_qual = quals->data + quals->at;
         quals->at += block_qual_len;
      }
      if (!keep) {
         for (uint64_t i = 0; i < ntags; ++i)
            str_next(rf, &r->tags, &rf->text);
         check(rf, rf->f);
         continue;
      }

      /* RNAME RNEXT PNEXT TLEN, as RAMFile and RAMNTupleRecord give them back. */
      for (size_t i = 0; i < ncigar; ++i) {
         if ((rf->cigar[i] & 0xf) > 8)
            die(rf, "bad CIGAR operation");
      }
      const int32_t tid = tid_of(rf, rf->tids, &rf->rnames, refid);
      int32_t mtid = tid_of(rf, rf->mtids, &rf->rnexts, refnext);
      if (mtid == kSame)
         mtid = tid;
      const int same_ref = tid >= 0 && mtid == tid;
      const int32_t mpos = same_ref ? wrap((int64_t)pnext_delta + pos) : pnext_delta;
      const int32_t tlen = wrap((int64_t)tlen_delta + predict_tlen(flag, pos, span, mpos, same_ref));
      if (flags & kSeqConsensus)
         decode_seq(rf, name, pos, (size_t)ncigar, (size_t)ndiffs);
      const Str *seq = flags & kSeqConsensus ? &rf->seq : &rf->seq_lit;
      const size_t l_seq = seq->n == 1 && seq->s[0] == '*' ? 0 : seq->n;

      /* QUAL without its offset: from the block, or from the record as RAMNTupleUtils::DecodeQuality()
         reads it. NULL is "*". */
      const char *qual = NULL;
      size_t qual_len = 0;
      if (flags & kQualInBlock) {
         qual = block_qual;
         qual_len = block_qual_len;
      } else if (flags & kDrop) {
         /* Dropped: none. */
      } else if (flags & kIlluminaBinning) {
         qual = rf->qual.n > 0 ? rf->qual.s : NULL;
         qual_len = rf->qual.n;
      } else if (rf->qual.n > 0 && !(rf->qual.n == 1 && rf->qual.s[0] == '*')) {
         for (size_t i = 0; i < rf->qual.n; ++i)
            rf->qual.s[i] = (char)(rf->qual.s[i] - '!');
         qual = rf->qual.s;
         qual_len = rf->qual.n;
      }
      if (qual != NULL && qual_len != l_seq)
         die(rf, "record %s has %zu bases and %zu qualities", name, l_seq, qual_len);
      if (bam_set1(b, name_len, name, flag, tid, pos, (uint8_t)mapq, (size_t)ncigar, rf->cigar, mtid, mpos, tlen, l_seq,
                   seq->s, qual, 0) < 0)
         die(rf, "record %s does not fit a bam1_t", name);

      /* The tags in input order, as RAMTagReader::Append() writes them. */
      uint64_t text_used = 0;
      const char *order = tag_order >= 0 && (uint64_t)tag_order < rf->orders.n ? rf->orders.v[tag_order] : "";
      const size_t order_len = strlen(order);
      for (size_t at = 0; at + 1 < order_len; at += 2) {
         int t = -1;
         for (int k = 0; k < 10 && t < 0; ++k)
            t = memcmp(order + at, kTypedTags[k], 2) == 0 ? k : -1;
         if (t < 0) {
            if (text_used < ntags) {
               str_next(rf, &r->tags, &rf->text);
               text_used++;
               text_tag(rf, b, rf->text.s, rf->text.n);
            }
            continue;
         }
         const char *tag = kTypedTags[t];
         if (t < kRG) {
            int_tag(rf, b, tag, tag_int[t]);
         } else if (t >= kMD) {
            aux(rf, b, tag, 'Z', rf->tag_str[t - kMD].n + 1, rf->tag_str[t - kMD].s);
         } else {
            const Strings *ids = t == kRG ? &rf->rg : &rf->pg;
            const int64_t id = t == kRG ? tag_rg : tag_pg;
            const char *v = id >= 0 && (uint64_t)id < ids->n ? ids->v[id] : "";
            aux(rf, b, tag, 'Z', strlen(v) + 1, v);
         }
      }
      /* Text tags the order does not account for, as in a damaged file, still come out. */
      for (; text_used < ntags; ++text_used) {
         str_next(rf, &r->tags, &rf->text);
         text_tag(rf, b, rf->text.s, rf->text.n);
      }
      check(rf, rf->f);
      return 1;
   }
   if (!rf->stopped && (names->at != names->ndata || quals->next != quals->nlengths))
      die(rf, "the last block has more names or qualities than records");
   return 0;
}

/* The setjmp of each public function, in a function of its own so that no variable lives across it. */
static int guarded_open(RamFile *rf, const char *path)
{
   if (setjmp(rf->trap))
      return -1;
   open_file(rf, path);
   return 0;
}

static int guarded_start(RamFile *rf, int query, int32_t qref, int32_t qstart, int32_t qend)
{
   if (setjmp(rf->trap))
      return -1;
   start(rf, query, qref, qstart, qend);
   return 0;
}

RamFile *ram_open(const char *path, char *err, size_t errlen)
{
   RamFile *rf = calloc(1, sizeof *rf);
   if (rf == NULL) {
      snprintf(err, errlen, "out of memory");
      return NULL;
   }
   if (guarded_open(rf, path) < 0) {
      snprintf(err, errlen, "%s", rf->err);
      ram_close(rf);
      return NULL;
   }
   return rf;
}

void ram_close(RamFile *rf)
{
   if (rf == NULL)
      return;
   records_free(&rf->r);
   block_free(&rf->names);
   block_free(&rf->quals);
   block_free(&rf->cons);
   rnt_cursor_free(&rf->cons_pos);
   rnt_cursor_free(&rf->mc);
   str_free(&rf->ms);
   if (rf->meta)
      rnt_close(rf->meta);
   if (rf->f)
      rnt_close(rf->f);
   sam_hdr_destroy(rf->hdr);
   Strings *all[] = {&rf->lines, &rf->rnames, &rf->rnexts, &rf->orders, &rf->rg, &rf->pg};
   for (size_t i = 0; i < sizeof all / sizeof *all; ++i)
      strings_free(all[i]);
   Str *strs[] = {&rf->seq,        &rf->seq_lit,    &rf->qual,       &rf->text,
                  &rf->tag_str[0], &rf->tag_str[1], &rf->tag_str[2], &rf->tag_str[3]};
   for (size_t i = 0; i < sizeof strs / sizeof *strs; ++i)
      free(strs[i]->s);
   free(rf->tids);
   free(rf->mtids);
   free(rf->ends);
   free(rf->cigar);
   free(rf->diffs);
   free(rf->array);
   free(rf->values);
   free(rf->varints);
   free(rf);
}

sam_hdr_t *ram_header(RamFile *rf)
{
   return rf->hdr;
}

int ram_query(RamFile *rf, const char *region)
{
   if (rf->failed)
      return -1;
   rf->err[0] = '\0';
   const int query = region != NULL && region[0] != '\0' && strcmp(region, "*") != 0;
   int32_t qref = -1;
   int32_t qstart = 0;
   int32_t qend = 0;
   if (query) {
      qref = parse_region(rf, region, &qstart, &qend);
      if (qref < 0)
         return -1;
   }
   return guarded_start(rf, query, qref, qstart, qend);
}

int ram_next(RamFile *rf, bam1_t *b)
{
   if (rf->failed)
      return -1;
   rf->err[0] = '\0';
   if (setjmp(rf->trap))
      return -1;
   return next_record(rf, b);
}

const char *ram_error(const RamFile *rf)
{
   return rf->err[0] ? rf->err : NULL;
}
