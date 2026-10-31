/* ramread: write a RAM file out as SAM without ROOT, as ramdump does.

   The record layout, the tag orders, the blocks of names and qualities and the
   PNEXT and TLEN predictions follow ramtools' RAMNTupleRecord, RAMTags and
   RecordBlocks; a change to those needs the same change here.

   Usage: ramread [-h] <file.ram> [region] */

#include "rntuple.h"

#include <htscodecs/fqzcomp_qual.h>
#include <htscodecs/rANS_static4x16.h>
#include <htscodecs/tokenise_name3.h>

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define die rnt_die

enum {
   kFormatVersion = 8
};

/* RAMNTupleRecord's compression_flags bits. */
enum {
   kIlluminaBinning = 1 << 15,
   kDrop = 1 << 16,
   kQualInBlock = 1 << 18
};

/* The first byte of a quality block: rANS or fqzcomp. */
enum {
   kQualRans = 1,
   kQualFqzcomp = 2
};

static void *xrealloc(void *p, size_t n)
{
   void *q = realloc(p, n ? n : 1);
   if (q == NULL)
      die("out of memory");
   return q;
}

/* A growable byte string. */
typedef struct {
   char *s;
   size_t n, cap;
} Str;

static void reserve(Str *b, size_t extra)
{
   if (extra > SIZE_MAX / 2 - b->n)
      die("out of memory");
   if (b->n + extra > b->cap) {
      b->cap = 2 * (b->n + extra);
      b->s = xrealloc(b->s, b->cap);
   }
}

static void put(Str *b, const void *s, size_t n)
{
   reserve(b, n);
   memcpy(b->s + b->n, s, n);
   b->n += n;
}

static void put_str(Str *b, const char *s)
{
   put(b, s, strlen(s));
}

static void put_char(Str *b, char c)
{
   put(b, &c, 1);
}

/* Decimal, without printf: a long-read CIGAR has thousands of numbers. */
static void put_int(Str *b, int64_t v)
{
   char tmp[24];
   size_t i = sizeof tmp;
   uint64_t u = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
   do {
      tmp[--i] = (char)('0' + u % 10);
      u /= 10;
   } while (u > 0);
   if (v < 0)
      tmp[--i] = '-';
   put(b, tmp + i, sizeof tmp - i);
}

/* A string field: its index column and its characters. */
typedef struct {
   RntCursor index, chars;
} StrField;

static void str_field(StrField *f, const RntFile *file, const char *path)
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

/* Reads the next string into b, replacing what it held. */
static void str_next(StrField *f, Str *b)
{
   const uint64_t n = rnt_len(&f->index);
   if (n > (1U << 30))
      die("string too long");
   b->n = 0;
   reserve(b, (size_t)n + 1);
   rnt_bytes(&f->chars, (size_t)n, (uint8_t *)b->s);
   b->n = (size_t)n;
   b->s[b->n] = '\0';
}

/* A std::vector<std::string> of METADATA's one entry. */
typedef struct {
   char **v;
   size_t n;
} Strings;

static Strings meta_strings(const RntFile *meta, const char *name)
{
   char path[64];
   snprintf(path, sizeof path, "%s._0", name);
   RntCursor index;
   rnt_cursor(&index, meta, rnt_column(meta, name, 0));
   StrField item;
   str_field(&item, meta, path);
   Strings out = {NULL, (size_t)rnt_len(&index)};
   if (out.n > (1U << 28))
      die("bad %s", name);
   out.v = xrealloc(NULL, out.n * sizeof *out.v);
   Str b = {0};
   for (size_t i = 0; i < out.n; ++i) {
      str_next(&item, &b);
      out.v[i] = xrealloc(NULL, b.n + 1);
      memcpy(out.v[i], b.s, b.n + 1);
   }
   free(b.s);
   rnt_cursor_free(&index);
   str_free(&item);
   return out;
}

static void strings_free(Strings *s)
{
   for (size_t i = 0; i < s->n; ++i)
      free(s->v[i]);
   free(s->v);
}

/* RAMFile::Name(): the name with that id, "*" for -1 or an id the file does not have. */
static const char *name_of(const Strings *s, int64_t id)
{
   return id >= 0 && (uint64_t)id < s->n ? s->v[id] : "*";
}

/* The typed tags of RAMTags.h, in the order of the record's tag_ members. */
static const char *const kTypedTags[10] = {"NM", "AS", "XS", "MQ", "RG", "PG", "MD", "MC", "SA", "XA"};
enum {
   kRG = 4,
   kMD = 6
};

typedef struct {
   RntCursor flag, refid, pos, mapq, cigar_n, cigar, refnext, pnext_delta, tlen_delta;
   StrField seq, qual;
   RntCursor tag_int[4], tag_rg, tag_pg;
   StrField tag_str[4];
   RntCursor tag_order, tags_n, flags;
   StrField tags;
} Records;

static void records_open(Records *r, const RntFile *f)
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
   str_field(&r->seq, f, "record.seq");
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
   StrField *strs[] = {&r->seq, &r->qual, &r->tag_str[0], &r->tag_str[1], &r->tag_str[2], &r->tag_str[3]};
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
                           &r->tag_order,  &r->tags_n,     &r->flags};
   for (size_t i = 0; i < sizeof cursors / sizeof *cursors; ++i)
      rnt_cursor_free(cursors[i]);
   StrField *strs[] = {&r->seq, &r->qual, &r->tag_str[0], &r->tag_str[1], &r->tag_str[2], &r->tag_str[3], &r->tags};
   for (size_t i = 0; i < sizeof strs / sizeof *strs; ++i)
      str_free(strs[i]);
}

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

/* Reads the blob of a block of rows, skipping the empty ones before it. */
static size_t read_blob(Block *b, uint64_t rows)
{
   size_t n = 0;
   for (uint64_t i = 0; i < rows; ++i) {
      const uint64_t len = rnt_len(&b->blob_n);
      if (len > (1U << 30))
         die("block too large");
      if (i + 1 < rows) {
         rnt_bytes(&b->blob, (size_t)len, NULL);
         continue;
      }
      n = (size_t)len;
      b->bytes = xrealloc(b->bytes, n);
      rnt_bytes(&b->blob, n, b->bytes);
   }
   return n;
}

/* NameBlockReader: the block's names, each ended by '\0'. */
static void load_names(Block *b, uint64_t rows)
{
   if (b->data && b->at != b->ndata)
      die("a block has more names than records");
   const size_t n = read_blob(b, rows);
   uint32_t size = 0;
   free(b->data);
   b->data = (char *)tok3_decode_names(b->bytes, (uint32_t)n, &size);
   if (b->data == NULL)
      die("the name tokeniser could not decode a block");
   uint64_t names = 0;
   for (uint32_t i = 0; i < size; ++i)
      names += b->data[i] == '\0';
   if (names != rows || (size > 0 && b->data[size - 1] != '\0'))
      die("a block has %" PRIu64 " names for %" PRIu64 " records", names, rows);
   b->ndata = size;
   b->at = 0;
}

static uint8_t *unrans(uint8_t *in, size_t n, unsigned *out)
{
   uint8_t *p = rans_uncompress_4x16(in, (unsigned)n, out);
   if (p == NULL)
      die("rANS could not decode a quality block");
   return p;
}

/* QualityBlockReader: the lossless qualities of the block's records that have kQualInBlock. */
static void load_quals(Block *b, uint64_t rows)
{
   if (b->next != b->nlengths)
      die("a block has more qualities than records");
   size_t records = 0;
   for (uint64_t i = 0; i < rows; ++i)
      records += (rnt_int(&b->flags) & kQualInBlock) != 0;
   const size_t n = read_blob(b, rows);
   free(b->data);
   b->data = NULL;
   b->ndata = b->at = b->next = 0;
   b->nlengths = records;
   b->lengths = xrealloc(b->lengths, records * sizeof *b->lengths);
   if (records == 0)
      return;
   if (n < 2)
      die("a quality block is damaged");
   uint8_t *in = b->bytes;
   size_t total = 0;
   if (in[0] == kQualRans) {
      /* The size of the lengths part, the lengths as varints, then the qualities: each part rANS. */
      if (n < 5)
         die("a quality block is damaged");
      const size_t nlen = (size_t)in[1] | (size_t)in[2] << 8 | (size_t)in[3] << 16 | (size_t)in[4] << 24;
      if (nlen >= n - 5)
         die("a quality block is damaged");
      unsigned nv = 0;
      uint8_t *v = unrans(in + 5, nlen, &nv);
      size_t k = 0;
      uint32_t x = 0;
      int shift = 0;
      for (unsigned i = 0; i < nv; ++i) {
         if (shift > 28)
            die("a quality block is damaged");
         x |= (uint32_t)(v[i] & 0x7f) << shift;
         shift += 7;
         if (v[i] < 0x80) {
            if (k == records || x > (uint32_t)INT32_MAX)
               die("a quality block is damaged");
            b->lengths[k++] = (int)x;
            total += x;
            x = 0;
            shift = 0;
         }
      }
      free(v);
      unsigned nq = 0;
      b->data = (char *)unrans(in + 5 + nlen, n - 5 - nlen, &nq);
      b->ndata = nq;
      if (shift != 0 || k != records || nq != total)
         die("a quality block is damaged");
   } else if (in[0] == kQualFqzcomp) {
      size_t size = 0;
      b->data = fqz_decompress((char *)in + 1, n - 1, &size, b->lengths, (int)records);
      if (b->data == NULL)
         die("fqzcomp could not decode a quality block");
      for (size_t i = 0; i < records; ++i) {
         if (b->lengths[i] < 0)
            die("a quality block is damaged");
         total += (size_t)b->lengths[i];
      }
      if (total != size)
         die("a quality block is damaged");
      b->ndata = size;
   } else {
      die("a quality block has an unknown codec");
   }
   for (size_t i = 0; i < b->ndata; ++i)
      b->data[i] = (char)(b->data[i] + '!');
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
   1-based and inclusive. Start and end come back 0-based. */
static int32_t parse_region(const char *region, const Strings *rnames, int32_t *start, int32_t *end)
{
   *start = 0;
   *end = INT32_MAX - 1;
   int32_t refid = find_ref(rnames, region, strlen(region));
   if (refid >= 0)
      return refid;
   const char *colon = strrchr(region, ':');
   const size_t name_len = colon ? (size_t)(colon - region) : strlen(region);
   if (colon) {
      const char *dash = strchr(colon + 1, '-');
      const int64_t a = position(colon + 1, dash ? (size_t)(dash - colon - 1) : strlen(colon + 1));
      const int64_t b = dash ? position(dash + 1, strlen(dash + 1)) : a;
      if (a < 0 || b < 0)
         die("invalid region %s: use rname[:start[-end]]", region);
      *start = a > 0 ? (int32_t)(a - 1) : 0;
      *end = (int32_t)(b - 1);
   }
   refid = find_ref(rnames, region, name_len);
   if (refid < 0)
      die("reference %.*s not found", (int)name_len, region);
   return refid;
}

/* The first row at or after (refid, pos) in a sorted file, where unplaced records sort last. */
static uint64_t first_row_at_or_after(const RntFile *f, uint64_t total, int32_t refid, int32_t pos)
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

int main(int argc, char **argv)
{
   int header = 0;
   int usage = argc < 2;
   const char *path = NULL;
   const char *region = NULL;
   for (int i = 1; i < argc; ++i) {
      if (strcmp(argv[i], "-h") == 0)
         header = 1;
      else if (path == NULL && argv[i][0] != '-')
         path = argv[i];
      else if (region == NULL && argv[i][0] != '-')
         region = argv[i];
      else
         usage = 1;
   }
   if (usage || path == NULL) {
      fputs("Usage: ramread [-h] <file.ram> [region]\n"
            "Region is rname[:start[-end]], 1-based and inclusive, as in samtools.\n",
            stderr);
      return 2;
   }

   RntFile *meta = rnt_open(path, "METADATA");
   RntCursor c;
   rnt_cursor(&c, meta, rnt_column(meta, "format_version", 0));
   const int64_t version = rnt_entries(meta) ? rnt_int(&c) : -1;
   rnt_cursor_free(&c);
   if (version != kFormatVersion)
      die("%s has RAM format version %" PRId64 "; this reader reads version %d", path, version, kFormatVersion);
   Strings header_lines = meta_strings(meta, "header_lines");
   Strings rnames = meta_strings(meta, "rname_refs");
   Strings rnexts = meta_strings(meta, "rnext_refs");
   Strings orders = meta_strings(meta, "tag_orders");
   Strings rg = meta_strings(meta, "tag_rg");
   Strings pg = meta_strings(meta, "tag_pg");
   rnt_cursor(&c, meta, rnt_column(meta, "max_ref_span", 0));
   const int32_t max_span = wrap(rnt_int(&c));
   rnt_cursor_free(&c);
   rnt_cursor(&c, meta, rnt_column(meta, "coordinate_sorted", 0));
   const int sorted = rnt_int(&c) != 0;
   rnt_cursor_free(&c);
   RntCursor ends_n;
   rnt_cursor(&ends_n, meta, rnt_column(meta, "block_ends", 0));
   rnt_cursor(&c, meta, rnt_column(meta, "block_ends._0", 0));
   const uint64_t nblocks = rnt_len(&ends_n);
   if (nblocks > (1U << 28))
      die("bad block_ends");
   uint64_t *ends = xrealloc(NULL, (size_t)nblocks * sizeof *ends);
   for (uint64_t i = 0; i < nblocks; ++i)
      ends[i] = (uint64_t)rnt_int(&c);
   rnt_cursor_free(&ends_n);
   rnt_cursor_free(&c);
   rnt_close(meta);

   RntFile *f = rnt_open(path, "RAM");
   const uint64_t entries = rnt_entries(f);
   for (uint64_t i = 0; i < nblocks; ++i) {
      if ((i > 0 && ends[i] <= ends[i - 1]) || ends[i] >= entries)
         die("the block ends in the file are damaged");
   }
   if (entries > 0 && (nblocks == 0 || ends[nblocks - 1] != entries - 1))
      die("records after the last block");

   /* No region, or "*", means every record. In a sorted file the scan starts a read's longest span before
      the region; otherwise it reads every record. */
   const int query = region != NULL && region[0] != '\0' && strcmp(region, "*") != 0;
   int32_t qref = -1;
   int32_t qstart = 0;
   int32_t qend = 0;
   uint64_t first = 0;
   if (query) {
      qref = parse_region(region, &rnames, &qstart, &qend);
      if (sorted)
         first = first_row_at_or_after(f, entries, qref, max_span > 0 && qstart > max_span ? qstart - max_span : 0);
   }
   /* The names and qualities come a block at a time, so decoding starts at the block that holds the first row. */
   uint64_t block = 0;
   while (block < nblocks && ends[block] < first)
      block++;
   const uint64_t block_start = block == 0 ? 0 : ends[block - 1] + 1;

   Records r;
   records_open(&r, f);
   Block names = {0};
   Block quals = {0};
   rnt_cursor(&names.blob_n, f, rnt_column(f, "nameblock", 0));
   rnt_cursor(&names.blob, f, rnt_column(f, "nameblock._0", 0));
   rnt_cursor(&quals.blob_n, f, rnt_column(f, "qualblock", 0));
   rnt_cursor(&quals.blob, f, rnt_column(f, "qualblock._0", 0));
   rnt_cursor(&quals.flags, f, rnt_column(f, "record.compression_flags", 0));
   if (block_start > 0 && block_start < entries) {
      records_seek(&r, block_start);
      rnt_seek(&names.blob_n, block_start);
      rnt_seek_items(&names.blob, &names.blob_n);
      rnt_seek(&quals.blob_n, block_start);
      rnt_seek_items(&quals.blob, &quals.blob_n);
      rnt_seek(&quals.flags, block_start);
   }

   Str line = {0};
   if (header) {
      for (size_t i = 0; i < header_lines.n; ++i) {
         put_str(&line, header_lines.v[i]);
         put_char(&line, '\n');
      }
   }

   Str seq = {0};
   Str qual = {0};
   Str text = {0};
   Str tag_str[4] = {{0}};
   uint32_t *cigar = NULL;
   size_t cigar_cap = 0;
   int stopped = 0;
   for (uint64_t row = block_start; row < entries && !stopped; ++row) {
      if (row == block_start || row == ends[block - 1] + 1) {
         const uint64_t rows = ends[block] - (block == 0 ? 0 : ends[block - 1] + 1) + 1;
         load_names(&names, rows);
         load_quals(&quals, rows);
         block++;
      }
      const uint16_t flag = (uint16_t)rnt_int(&r.flag);
      const int64_t refid = rnt_int(&r.refid);
      const int32_t pos = (int32_t)rnt_int(&r.pos);
      const int64_t mapq = rnt_int(&r.mapq);
      const uint64_t ncigar = rnt_len(&r.cigar_n);
      if (ncigar > (1U << 28))
         die("CIGAR too long");
      if (ncigar > cigar_cap) {
         cigar_cap = (size_t)ncigar;
         cigar = xrealloc(cigar, cigar_cap * sizeof *cigar);
      }
      for (size_t i = 0; i < ncigar; ++i)
         cigar[i] = (uint32_t)rnt_int(&r.cigar);
      const int64_t refnext = rnt_int(&r.refnext);
      const int32_t pnext_delta = (int32_t)rnt_int(&r.pnext_delta);
      const int32_t tlen_delta = (int32_t)rnt_int(&r.tlen_delta);
      str_next(&r.seq, &seq);
      str_next(&r.qual, &qual);
      int32_t tag_int[4];
      for (int i = 0; i < 4; ++i)
         tag_int[i] = (int32_t)rnt_int(&r.tag_int[i]);
      const int64_t tag_rg = rnt_int(&r.tag_rg);
      const int64_t tag_pg = rnt_int(&r.tag_pg);
      for (int i = 0; i < 4; ++i)
         str_next(&r.tag_str[i], &tag_str[i]);
      const int64_t tag_order = rnt_int(&r.tag_order);
      const uint64_t ntags = rnt_len(&r.tags_n);
      const uint32_t flags = (uint32_t)rnt_int(&r.flags);

      /* Whether ramntuplescan() keeps the row. The rest of the row is decoded either way: the cursors and the
         blocks move one row at a time. */
      int keep = 1;
      if (query) {
         /* A record with no reference span still covers the base it sits on, as htslib's bam_endpos says. */
         const int64_t span = ref_span(cigar, (size_t)ncigar);
         keep = refid == qref && pos <= qend && (pos >= qstart || pos + (span > 0 ? span : 1) - 1 >= qstart);
         if (sorted && row >= first && (refid < 0 || refid > qref || (refid == qref && pos > qend)))
            stopped = 1;
         keep = keep && row >= first && !stopped;
      }
      const size_t line_start = line.n;

      /* QNAME, from the block. */
      const char *name = names.data + names.at;
      const size_t name_len = strlen(name);
      names.at += name_len + 1;
      put(&line, name, name_len);

      /* FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN, as RAMFile and RAMNTupleRecord give them back. */
      const char *rname = name_of(&rnames, refid);
      const char *rnext = name_of(&rnexts, refnext);
      const int same_ref = strcmp(rname, "*") != 0 && (strcmp(rnext, "=") == 0 || strcmp(rnext, rname) == 0);
      const int32_t pnext = (same_ref ? wrap((int64_t)pnext_delta + pos) : pnext_delta) + 1;
      const int64_t span = ref_span(cigar, (size_t)ncigar);
      const int32_t tlen = wrap((int64_t)tlen_delta + predict_tlen(flag, pos, span, pnext - 1, same_ref));
      put_char(&line, '\t');
      put_int(&line, flag);
      put_char(&line, '\t');
      put_str(&line, rname);
      put_char(&line, '\t');
      put_int(&line, (int64_t)pos + 1);
      put_char(&line, '\t');
      put_int(&line, mapq);
      put_char(&line, '\t');
      if (ncigar == 0)
         put_char(&line, '*');
      for (size_t i = 0; i < ncigar; ++i) {
         if ((cigar[i] & 0xf) > 8)
            die("bad CIGAR operation");
         put_int(&line, cigar[i] >> 4);
         put_char(&line, "MIDNSHP=X"[cigar[i] & 0xf]);
      }
      put_char(&line, '\t');
      put_str(&line, rnext);
      put_char(&line, '\t');
      put_int(&line, pnext);
      put_char(&line, '\t');
      put_int(&line, tlen);
      put_char(&line, '\t');
      put(&line, seq.s, seq.n);
      put_char(&line, '\t');

      /* QUAL: from the block, or from the record as RAMNTupleUtils::DecodeQuality() reads it. */
      if (flags & kQualInBlock) {
         if (quals.next == quals.nlengths)
            die("more qualities in records than in their block");
         const size_t len = (size_t)quals.lengths[quals.next++];
         put(&line, quals.data + quals.at, len);
         quals.at += len;
      } else if (flags & kDrop) {
         put_char(&line, '*');
      } else if (flags & kIlluminaBinning) {
         if (qual.n == 0)
            put_char(&line, '*');
         for (size_t i = 0; i < qual.n; ++i)
            put_char(&line, (char)(qual.s[i] + 33));
      } else {
         put(&line, qual.s, qual.n);
      }

      /* The tags in input order, as RAMTagReader::Append() writes them. */
      uint64_t text_used = 0;
      const char *order = tag_order >= 0 && (uint64_t)tag_order < orders.n ? orders.v[tag_order] : "";
      const size_t order_len = strlen(order);
      for (size_t at = 0; at + 1 < order_len; at += 2) {
         int t = -1;
         for (int k = 0; k < 10 && t < 0; ++k)
            t = memcmp(order + at, kTypedTags[k], 2) == 0 ? k : -1;
         if (t < 0) {
            if (text_used < ntags) {
               str_next(&r.tags, &text);
               text_used++;
               put_char(&line, '\t');
               put(&line, text.s, text.n);
            }
            continue;
         }
         put_char(&line, '\t');
         put_str(&line, kTypedTags[t]);
         put_str(&line, t < kRG ? ":i:" : ":Z:");
         if (t < kRG)
            put_int(&line, tag_int[t]);
         else if (t >= kMD)
            put(&line, tag_str[t - kMD].s, tag_str[t - kMD].n);
         else if (t == kRG)
            put_str(&line, tag_rg >= 0 && (uint64_t)tag_rg < rg.n ? rg.v[tag_rg] : "");
         else
            put_str(&line, tag_pg >= 0 && (uint64_t)tag_pg < pg.n ? pg.v[tag_pg] : "");
      }
      /* Text tags the order does not account for, as in a damaged file, still come out. */
      for (; text_used < ntags; ++text_used) {
         str_next(&r.tags, &text);
         put_char(&line, '\t');
         put(&line, text.s, text.n);
      }
      put_char(&line, '\n');
      if (!keep)
         line.n = line_start;
      if (line.n >= (1U << 20)) {
         fwrite(line.s, 1, line.n, stdout);
         line.n = 0;
      }
   }
   if (line.n > 0)
      fwrite(line.s, 1, line.n, stdout);
   if (fflush(stdout) != 0)
      die("cannot write the output");
   if (!stopped && (names.at != names.ndata || quals.next != quals.nlengths))
      die("the last block has more names or qualities than records");

   free(line.s);
   free(seq.s);
   free(qual.s);
   free(text.s);
   for (int i = 0; i < 4; ++i)
      free(tag_str[i].s);
   free(cigar);
   free(ends);
   records_free(&r);
   block_free(&names);
   block_free(&quals);
   rnt_close(f);
   Strings *all[] = {&header_lines, &rnames, &rnexts, &orders, &rg, &pg};
   for (size_t i = 0; i < sizeof all / sizeof *all; ++i)
      strings_free(all[i]);
   return 0;
}
