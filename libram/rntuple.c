#include "rntuple.h"

#define XXH_INLINE_ALL
#include <xxhash.h>
#include <zstd.h>

#include <inttypes.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static _Noreturn void die(RntFile *f, const char *fmt, ...);

static void *xrealloc(RntFile *f, void *p, size_t n)
{
   void *q = realloc(p, n ? n : 1);
   if (q == NULL)
      die(f, "out of memory");
   return q;
}

/* A byte range being parsed. Every read is checked against its end: the file is untrusted. */
typedef struct {
   const uint8_t *p;
   size_t n, off;
   RntFile *f;
} Buf;

static const uint8_t *take(Buf *c, size_t n)
{
   if (n > c->n - c->off)
      die(c->f, "truncated record");
   const uint8_t *q = c->p + c->off;
   c->off += n;
   return q;
}

static void seek(Buf *c, uint64_t off)
{
   if (off > c->n)
      die(c->f, "offset past the end of a record");
   c->off = (size_t)off;
}

static uint64_t le(Buf *c, size_t n)
{
   const uint8_t *q = take(c, n);
   uint64_t v = 0;
   for (size_t i = n; i-- > 0;)
      v = (v << 8) | q[i];
   return v;
}

static uint64_t be(Buf *c, size_t n)
{
   const uint8_t *q = take(c, n);
   uint64_t v = 0;
   for (size_t i = 0; i < n; ++i)
      v = (v << 8) | q[i];
   return v;
}

static uint64_t le64_at(const uint8_t *q)
{
   uint64_t v = 0;
   for (size_t i = 8; i-- > 0;)
      v = (v << 8) | q[i];
   return v;
}

/* A ROOT TString: one length byte, or 255 followed by a 32-bit big-endian length. */
static void tstring(Buf *c, const uint8_t **s, size_t *len)
{
   size_t n = (size_t)be(c, 1);
   if (n == 255)
      n = (size_t)be(c, 4);
   *s = take(c, n);
   *len = n;
}

/* An RNTuple string: a 32-bit little-endian length and the bytes. */
static char *ntstring(Buf *c)
{
   const size_t n = (size_t)le(c, 4);
   const uint8_t *s = take(c, n);
   char *out = xrealloc(c->f, NULL, n + 1);
   memcpy(out, s, n);
   out[n] = '\0';
   return out;
}

/* Opens a frame and returns the offset where it ends; a list frame also gives its number of items. */
static size_t frame(Buf *c, int want_list, uint32_t *nitems)
{
   const size_t start = c->off;
   const int64_t size = (int64_t)le(c, 8);
   const int is_list = size < 0;
   const uint64_t len = is_list ? (uint64_t)0 - (uint64_t)size : (uint64_t)size;
   if (is_list != want_list)
      die(c->f, "expected a %s frame", want_list ? "list" : "record");
   if (is_list)
      *nitems = (uint32_t)le(c, 4);
   if (len > c->n - start || start + len < c->off)
      die(c->f, "bad frame size");
   return start + (size_t)len;
}

typedef struct {
   uint64_t offset, size; /* the locator, without the checksum after it */
   uint64_t first;        /* index in the column of the page's first element */
   uint32_t nelem, cluster;
   int checksum;
} Page;

struct RntColumn {
   uint16_t type;
   uint16_t repr;
   uint32_t field;
   Page *pages;
   size_t npages;
   uint64_t *cluster_first; /* index in the column of each cluster's first element */
   uint32_t nclusters;
};

struct RntFile {
   FILE *f;
   uint64_t size;
   jmp_buf *trap; /* where an error unwinds to */
   char err[256]; /* the first error, or empty */
   void *tmp[3];  /* buffers that live only while an envelope or a page is read */
   char **names;
   uint32_t *parents;
   uint32_t nfields;
   RntColumn *columns;
   uint32_t ncolumns;
   uint64_t entries;
};

/* Errors unwind with longjmp to the trap that rnt_open or the loading of a page set. Every buffer belongs to
   the file, to a cursor or to the file's temporaries, so nothing leaks on the way. */
static void set_error(RntFile *f, const char *fmt, va_list ap)
{
   if (f->err[0] == '\0')
      vsnprintf(f->err, sizeof f->err, fmt, ap);
}

static _Noreturn void die(RntFile *f, const char *fmt, ...)
{
   va_list ap;
   va_start(ap, fmt);
   set_error(f, fmt, ap);
   va_end(ap);
   if (f->trap == NULL)
      abort();
   longjmp(*f->trap, 1);
}

/* Records the first error without unwinding; the cursors return zeros from then on. */
static void fail(RntFile *f, const char *fmt, ...)
{
   va_list ap;
   va_start(ap, fmt);
   set_error(f, fmt, ap);
   va_end(ap);
}

static void *temporary(RntFile *f, size_t n)
{
   for (size_t i = 0; i < sizeof f->tmp / sizeof *f->tmp; ++i) {
      if (f->tmp[i] == NULL)
         return f->tmp[i] = xrealloc(f, NULL, n);
   }
   abort();
}

/* Frees a temporary, or with keep, hands it to the caller. */
static void release(RntFile *f, void *p, int keep)
{
   for (size_t i = 0; i < sizeof f->tmp / sizeof *f->tmp; ++i) {
      if (f->tmp[i] == p && p != NULL)
         f->tmp[i] = NULL;
   }
   if (!keep)
      free(p);
}

static void free_temporaries(RntFile *f)
{
   for (size_t i = 0; i < sizeof f->tmp / sizeof *f->tmp; ++i) {
      free(f->tmp[i]);
      f->tmp[i] = NULL;
   }
}

static uint8_t *read_at(RntFile *f, uint64_t off, uint64_t n)
{
   if (off > f->size || n > f->size - off)
      die(f, "range %" PRIu64 "+%" PRIu64 " is past the end of the file", off, n);
   uint8_t *buf = temporary(f, (size_t)n);
   if (fseeko(f->f, (off_t)off, SEEK_SET) != 0 || fread(buf, 1, (size_t)n, f->f) != n)
      die(f, "cannot read %" PRIu64 " bytes at %" PRIu64, n, off);
   return buf;
}

/* Undoes ROOT's compression blocks: chunks of a 9-byte header and the compressed bytes. */
static uint8_t *unzip(RntFile *f, uint8_t *in, size_t nin, size_t nout)
{
   if (nin == nout)
      return in;
   uint8_t *out = temporary(f, nout);
   size_t i = 0;
   size_t o = 0;
   while (o < nout) {
      if (nin - i < 9)
         die(f, "truncated compression block");
      const uint8_t *h = in + i;
      const size_t csize = (size_t)h[3] | (size_t)h[4] << 8 | (size_t)h[5] << 16;
      const size_t usize = (size_t)h[6] | (size_t)h[7] << 8 | (size_t)h[8] << 16;
      if (h[0] != 'Z' || h[1] != 'S')
         die(f, "compression %c%c is not supported, only zstd", h[0], h[1]);
      if (csize > nin - i - 9 || usize > nout - o)
         die(f, "bad compression block");
      if (ZSTD_decompress(out + o, usize, h + 9, csize) != usize)
         die(f, "zstd failed");
      i += 9 + csize;
      o += usize;
   }
   if (i != nin)
      die(f, "bad compression block");
   release(f, in, 0);
   return out;
}

/* Reads an envelope, checks its type and checksum, and leaves the cursor on its payload. */
static uint8_t *
envelope(RntFile *f, uint64_t at, uint64_t nbytes, uint64_t len, unsigned type, Buf *c, uint64_t *checksum)
{
   if (len < 16 || len > SIZE_MAX)
      die(f, "bad envelope length");
   uint8_t *buf = unzip(f, read_at(f, at, nbytes), (size_t)nbytes, (size_t)len);
   const uint64_t ck = le64_at(buf + len - 8);
   if (XXH3_64bits(buf, (size_t)len - 8) != ck)
      die(f, "envelope %u checksum mismatch", type);
   *c = (Buf){buf, (size_t)len - 8, 0, f};
   const uint64_t id = le(c, 8);
   if ((id & 0xffff) != type || id >> 16 != len)
      die(f, "expected envelope %u of %" PRIu64 " bytes", type, len);
   if (checksum)
      *checksum = ck;
   return buf;
}

typedef struct {
   uint64_t size, offset;
} Locator;

static Locator locator(Buf *c)
{
   const int32_t head = (int32_t)le(c, 4);
   if (head >= 0)
      return (Locator){(uint64_t)head, le(c, 8)};
   /* Non-standard: 16 bits of locator size, 8 reserved, then the type as a negative 8-bit number. */
   const int type = -(int)(int8_t)(uint8_t)((uint32_t)head >> 24);
   if (type != 1 || ((uint32_t)head & 0xffff) != 20)
      die(c->f, "locator type %d is not supported", type);
   const uint64_t size = le(c, 8);
   return (Locator){size, le(c, 8)};
}

static void check_features(Buf *c)
{
   if (le(c, 8) != 0)
      die(c->f, "the file uses RNTuple features this reader does not know");
}

/* The four list frames of a schema description, in the header or its extension in the footer. */
static void schema(Buf *c, RntFile *f)
{
   uint32_t n = 0;
   size_t end = frame(c, 1, &n);
   if (n > c->n)
      die(f, "bad field count");
   f->names = xrealloc(f, f->names, (f->nfields + n) * sizeof *f->names);
   f->parents = xrealloc(f, f->parents, (f->nfields + n) * sizeof *f->parents);
   for (uint32_t i = 0; i < n; ++i) {
      const size_t fend = frame(c, 0, NULL);
      take(c, 8); /* field and type version */
      f->parents[f->nfields] = (uint32_t)le(c, 4);
      take(c, 4); /* structural role, flags */
      f->names[f->nfields++] = ntstring(c);
      c->off = fend;
   }
   c->off = end;

   end = frame(c, 1, &n);
   if (n > c->n)
      die(f, "bad column count");
   f->columns = xrealloc(f, f->columns, (f->ncolumns + n) * sizeof *f->columns);
   for (uint32_t i = 0; i < n; ++i) {
      const size_t cend = frame(c, 0, NULL);
      RntColumn *col = &f->columns[f->ncolumns++];
      memset(col, 0, sizeof *col);
      col->type = (uint16_t)le(c, 2);
      take(c, 2); /* bits on storage */
      col->field = (uint32_t)le(c, 4);
      if (le(c, 2) & 1)
         die(f, "deferred columns are not supported");
      col->repr = (uint16_t)le(c, 2);
      c->off = cend;
   }
   c->off = end;

   for (int k = 0; k < 2; ++k) { /* alias columns, extra type information */
      end = frame(c, 1, &n);
      c->off = end;
   }
}

/* Finds the key of the RNTuple called name in the top directory and returns the offset of its object. */
static uint64_t find_anchor(RntFile *f, const char *name, uint64_t *objlen)
{
   uint8_t *hb = read_at(f, 0, f->size < 64 ? f->size : 64);
   Buf c = {hb, f->size < 64 ? (size_t)f->size : 64, 0, f};
   if (memcmp(take(&c, 4), "root", 4) != 0)
      die(f, "not a ROOT file");
   take(&c, 4); /* version */
   const uint64_t begin = be(&c, 4);
   release(f, hb, 0);
   if (begin >= f->size)
      die(f, "bad file header");

   /* The top directory: its key, the file's name and title, then the streamed TDirectory. */
   uint8_t *db = read_at(f, begin, f->size - begin < 512 ? f->size - begin : 512);
   c = (Buf){db, f->size - begin < 512 ? (size_t)(f->size - begin) : 512, 0, f};
   take(&c, 14);
   seek(&c, be(&c, 2)); /* key length */
   const uint8_t *s = NULL;
   size_t len = 0;
   tstring(&c, &s, &len);
   tstring(&c, &s, &len);
   const unsigned dir_version = (unsigned)be(&c, 2);
   take(&c, 8); /* creation and modification time */
   const uint64_t nbytes_keys = be(&c, 4);
   take(&c, 4); /* nbytes name */
   const size_t w = dir_version > 1000 ? 8 : 4;
   take(&c, 2 * w); /* seek dir, seek parent */
   const uint64_t seek_keys = be(&c, w);
   release(f, db, 0);

   /* The key list: a key of its own, the number of keys, then each key's header. */
   uint8_t *kb = read_at(f, seek_keys, nbytes_keys);
   c = (Buf){kb, (size_t)nbytes_keys, 0, f};
   take(&c, 14);
   seek(&c, be(&c, 2));
   const uint64_t nkeys = be(&c, 4);
   uint64_t at = 0;
   uint64_t best_cycle = 0;
   for (uint64_t i = 0; i < nkeys; ++i) {
      take(&c, 4); /* nbytes */
      const unsigned version = (unsigned)be(&c, 2);
      const uint64_t len_obj = be(&c, 4);
      take(&c, 4); /* date */
      take(&c, 2); /* key length */
      const uint64_t cycle = be(&c, 2);
      const uint64_t seek_key = be(&c, version > 1000 ? 8 : 4);
      take(&c, version > 1000 ? 8 : 4); /* seek parent directory */
      const uint8_t *cls = NULL;
      const uint8_t *key = NULL;
      size_t ncls = 0;
      size_t nkey = 0;
      tstring(&c, &cls, &ncls);
      tstring(&c, &key, &nkey);
      tstring(&c, &s, &len);
      if (ncls == 13 && memcmp(cls, "ROOT::RNTuple", 13) == 0 && nkey == strlen(name) && memcmp(key, name, nkey) == 0 &&
          cycle > best_cycle) {
         best_cycle = cycle;
         at = seek_key;
         *objlen = len_obj;
      }
   }
   release(f, kb, 0);
   if (at == 0)
      die(f, "no RNTuple called %s", name);

   /* The object follows the key's own header. */
   uint8_t *ab = read_at(f, at, 18);
   c = (Buf){ab, 18, 0, f};
   const uint64_t nbytes = be(&c, 4);
   take(&c, 10);
   const uint64_t keylen = be(&c, 2);
   release(f, ab, 0);
   if (nbytes < keylen || nbytes - keylen != *objlen)
      die(f, "compressed anchors are not supported");
   return at + keylen;
}

/* Adds the pages of every column of the clusters in one page list envelope. */
static void page_list(RntFile *f, Buf *p, uint64_t header_ck, uint32_t first_cluster)
{
   if (le(p, 8) != header_ck)
      die(f, "page list belongs to another header");
   uint32_t n = 0;
   size_t end = frame(p, 1, &n);
   for (uint32_t i = 0; i < n; ++i) {
      const size_t send = frame(p, 0, NULL);
      take(p, 8); /* first entry */
      const uint64_t entries = le(p, 8);
      if (entries >> 56 & 1)
         die(f, "sharded clusters are not supported");
      f->entries += entries & ((1ULL << 56) - 1);
      p->off = send;
   }
   p->off = end;

   uint32_t nclusters = 0;
   const size_t clusters_end = frame(p, 1, &nclusters);
   for (uint32_t k = 0; k < nclusters; ++k) {
      uint32_t ncols = 0;
      const size_t cols_end = frame(p, 1, &ncols);
      if (ncols != f->ncolumns)
         die(f, "a cluster has %u columns, the schema %u", ncols, f->ncolumns);
      for (uint32_t i = 0; i < ncols; ++i) {
         RntColumn *col = &f->columns[i];
         if (col->nclusters != first_cluster + k)
            die(f, "page lists out of order");
         uint32_t npages = 0;
         const size_t pages_end = frame(p, 1, &npages);
         if (npages > p->n)
            die(f, "bad page count");
         col->pages = xrealloc(f, col->pages, (col->npages + npages) * sizeof *col->pages);
         Page *first = col->pages + col->npages;
         uint64_t elements = 0;
         for (uint32_t pg = 0; pg < npages; ++pg) {
            const int32_t nelem = (int32_t)le(p, 4);
            const Locator loc = locator(p);
            Page *page = &col->pages[col->npages++];
            page->checksum = nelem < 0;
            page->nelem = nelem < 0 ? (uint32_t)0 - (uint32_t)nelem : (uint32_t)nelem;
            page->offset = loc.offset;
            page->size = loc.size;
            page->cluster = first_cluster + k;
            page->first = elements;
            elements += page->nelem;
         }
         const int64_t offset = (int64_t)le(p, 8);
         if (offset < 0)
            die(f, "suppressed columns are not supported");
         col->cluster_first = xrealloc(f, col->cluster_first, (col->nclusters + 1) * sizeof *col->cluster_first);
         col->cluster_first[col->nclusters++] = (uint64_t)offset;
         for (uint32_t pg = 0; pg < npages; ++pg)
            first[pg].first += (uint64_t)offset;
         p->off = pages_end;
      }
      p->off = cols_end;
   }
   p->off = clusters_end;
}

static void parse(RntFile *f, const char *path, const char *name)
{
   f->f = fopen(path, "rb");
   struct stat st;
   if (f->f == NULL || fstat(fileno(f->f), &st) != 0)
      die(f, "cannot open %s", path);
   f->size = (uint64_t)st.st_size;

   /* The anchor: big-endian, with an XXH3 over everything after byte count and class version. */
   uint64_t objlen = 0;
   const uint64_t at = find_anchor(f, name, &objlen);
   if (objlen < 78)
      die(f, "anchor too short");
   uint8_t *ab = read_at(f, at, objlen);
   Buf c = {ab, (size_t)objlen, 0, f};
   seek(&c, objlen - 8);
   if (XXH3_64bits(ab + 6, (size_t)objlen - 14) != be(&c, 8))
      die(f, "anchor checksum mismatch");
   seek(&c, 6);
   if (be(&c, 2) != 1)
      die(f, "RNTuple epoch is not 1");
   take(&c, 6); /* major, minor, patch */
   const uint64_t seek_header = be(&c, 8);
   const uint64_t nbytes_header = be(&c, 8);
   const uint64_t len_header = be(&c, 8);
   const uint64_t seek_footer = be(&c, 8);
   const uint64_t nbytes_footer = be(&c, 8);
   const uint64_t len_footer = be(&c, 8);
   release(f, ab, 0);

   uint64_t header_ck = 0;
   uint8_t *hb = envelope(f, seek_header, nbytes_header, len_header, 1, &c, &header_ck);
   check_features(&c);
   for (int i = 0; i < 3; ++i)
      free(ntstring(&c)); /* name, description, writer */
   schema(&c, f);
   release(f, hb, 0);

   uint8_t *fb = envelope(f, seek_footer, nbytes_footer, len_footer, 2, &c, NULL);
   check_features(&c);
   if (le(&c, 8) != header_ck)
      die(f, "footer belongs to another header");
   const size_t ext_end = frame(&c, 0, NULL);
   if (c.off < ext_end)
      schema(&c, f);
   c.off = ext_end;
   for (uint32_t i = 0; i < f->ncolumns; ++i) {
      if (f->columns[i].field >= f->nfields || f->columns[i].repr != 0)
         die(f, "column %u is not supported", i);
   }

   uint32_t ngroups = 0;
   uint32_t clusters = 0;
   const size_t groups_end = frame(&c, 1, &ngroups);
   for (uint32_t g = 0; g < ngroups; ++g) {
      const size_t gend = frame(&c, 0, NULL);
      take(&c, 16); /* min entry, entry span */
      const uint32_t nclusters = (uint32_t)le(&c, 4);
      const uint64_t len_pl = le(&c, 8);
      const Locator pl = locator(&c);
      c.off = gend;
      Buf p;
      uint8_t *pb = envelope(f, pl.offset, pl.size, len_pl, 3, &p, NULL);
      page_list(f, &p, header_ck, clusters);
      release(f, pb, 0);
      clusters += nclusters;
   }
   c.off = groups_end;
   release(f, fb, 0);
}

RntFile *rnt_open(const char *path, const char *name, char *err, size_t errlen)
{
   /* volatile: GCC cannot see that f stays the same across the longjmp. */
   RntFile *volatile f = calloc(1, sizeof *f);
   if (f == NULL) {
      snprintf(err, errlen, "out of memory");
      return NULL;
   }
   jmp_buf trap;
   f->trap = &trap;
   if (setjmp(trap) != 0) {
      snprintf(err, errlen, "%s", f->err);
      rnt_close(f);
      return NULL;
   }
   parse(f, path, name);
   f->trap = NULL;
   return f;
}

void rnt_close(RntFile *f)
{
   free_temporaries(f);
   for (uint32_t i = 0; i < f->nfields; ++i)
      free(f->names[i]);
   for (uint32_t i = 0; i < f->ncolumns; ++i) {
      free(f->columns[i].pages);
      free(f->columns[i].cluster_first);
   }
   free(f->names);
   free(f->parents);
   free(f->columns);
   if (f->f)
      fclose(f->f);
   free(f);
}

uint64_t rnt_entries(const RntFile *f)
{
   return f->entries;
}

const char *rnt_error(const RntFile *f)
{
   return f->err[0] ? f->err : NULL;
}

const RntColumn *rnt_column(RntFile *f, const char *path, unsigned which)
{
   int64_t parent = -1;
   for (const char *p = path; *p;) {
      const char *dot = strchr(p, '.');
      const size_t len = dot ? (size_t)(dot - p) : strlen(p);
      int64_t found = -1;
      for (uint32_t i = 0; i < f->nfields && found < 0; ++i) {
         const int top = f->parents[i] == i;
         const int under = parent < 0 ? top : (!top && f->parents[i] == (uint64_t)parent);
         if (under && strlen(f->names[i]) == len && memcmp(f->names[i], p, len) == 0)
            found = i;
      }
      if (found < 0) {
         fail(f, "no field %s", path);
         return NULL;
      }
      parent = found;
      p += len + (dot != NULL);
   }
   for (uint32_t i = 0; i < f->ncolumns; ++i) {
      if (f->columns[i].field == (uint64_t)parent && which-- == 0)
         return &f->columns[i];
   }
   fail(f, "field %s has no such column", path);
   return NULL;
}

/* How the column type stores an element: its size in memory, and whether it is split, zigzag, delta, signed,
   an index or a bit. */
typedef struct {
   size_t size;
   int split, zigzag, delta, is_signed, index, bit;
} Encoding;

static Encoding encoding(uint16_t type)
{
   switch (type) {
   case 0x00: return (Encoding){1, 0, 0, 0, 0, 0, 1};
   case 0x01:
   case 0x02:
   case 0x04: return (Encoding){1, 0, 0, 0, 0, 0, 0};
   case 0x03: return (Encoding){1, 0, 0, 0, 1, 0, 0};
   case 0x05: return (Encoding){2, 0, 0, 0, 1, 0, 0};
   case 0x06: return (Encoding){2, 0, 0, 0, 0, 0, 0};
   case 0x07: return (Encoding){4, 0, 0, 0, 1, 0, 0};
   case 0x08: return (Encoding){4, 0, 0, 0, 0, 0, 0};
   case 0x09: return (Encoding){8, 0, 0, 0, 1, 0, 0};
   case 0x0A: return (Encoding){8, 0, 0, 0, 0, 0, 0};
   case 0x0E: return (Encoding){4, 0, 0, 0, 0, 1, 0};
   case 0x0F: return (Encoding){8, 0, 0, 0, 0, 1, 0};
   case 0x11: return (Encoding){2, 1, 1, 0, 1, 0, 0};
   case 0x12: return (Encoding){2, 1, 0, 0, 0, 0, 0};
   case 0x13: return (Encoding){4, 1, 1, 0, 1, 0, 0};
   case 0x14: return (Encoding){4, 1, 0, 0, 0, 0, 0};
   case 0x15: return (Encoding){8, 1, 1, 0, 1, 0, 0};
   case 0x16: return (Encoding){8, 1, 0, 0, 0, 0, 0};
   case 0x1A: return (Encoding){4, 1, 0, 1, 0, 1, 0};
   case 0x1B: return (Encoding){8, 1, 0, 1, 0, 1, 0};
   default: return (Encoding){0, 0, 0, 0, 0, 0, 0};
   }
}

/* What a cursor reads when its column is missing or not supported: a column of bytes that has no pages. */
static const RntColumn kNoColumn = {0x01, 0, 0, NULL, 0, NULL, 0};

void rnt_cursor(RntCursor *c, RntFile *f, const RntColumn *col)
{
   memset(c, 0, sizeof *c);
   c->file = f;
   c->col = &kNoColumn;
   if (col == NULL)
      fail(f, "no such column");
   else if (encoding(col->type).size == 0)
      fail(f, "column type 0x%02x is not supported", col->type);
   else
      c->col = col;
}

void rnt_cursor_free(RntCursor *c)
{
   free(c->bytes);
   free(c->ints);
   memset(c, 0, sizeof *c);
}

/* Reads the next page of the cursor's column and decodes its elements. */
static void read_page(RntCursor *c)
{
   RntFile *f = c->file;
   const RntColumn *col = c->col;
   if (c->page == col->npages)
      die(f, "a column ends before its last entry");
   const Page *pg = &col->pages[c->page];
   const Encoding e = encoding(col->type);
   c->consumed += c->n;
   if (pg->first != c->consumed)
      die(f, "the pages of a column do not follow each other");
   if (c->page == 0 || pg->cluster != c->cluster)
      c->prev = 0;
   c->cluster = pg->cluster;
   c->page++;

   if (e.size == 0 || pg->nelem > (1U << 30) / e.size || pg->size > UINT64_MAX - 8)
      die(f, "bad page size");
   const size_t n = pg->nelem;
   uint8_t *raw = read_at(f, pg->offset, pg->size + (pg->checksum ? 8 : 0));
   if (pg->checksum && XXH3_64bits(raw, (size_t)pg->size) != le64_at(raw + pg->size))
      die(f, "page checksum mismatch at %" PRIu64, pg->offset);
   uint8_t *page = unzip(f, raw, (size_t)pg->size, e.bit ? (n + 7) / 8 : n * e.size);
   c->n = n;
   c->at = 0;
   if (e.bit) {
      /* Eight to a byte, the first in the lowest bit. */
      c->bytes = xrealloc(f, c->bytes, n);
      for (size_t i = 0; i < n; ++i)
         c->bytes[i] = (uint8_t)(page[i / 8] >> (i % 8) & 1);
      release(f, page, 0);
      return;
   }
   if (e.size == 1) {
      free(c->bytes);
      release(f, page, 1);
      c->bytes = page;
      return;
   }
   c->ints = xrealloc(f, c->ints, (n ? n : 1) * sizeof *c->ints);
   const unsigned bits = (unsigned)(8 * e.size);
   for (size_t i = 0; i < n; ++i) {
      uint64_t u = 0;
      for (size_t b = e.size; b-- > 0;)
         u = (u << 8) | page[e.split ? b * n + i : i * e.size + b];
      int64_t v = (int64_t)u;
      if (e.zigzag)
         v = (int64_t)(u >> 1) ^ -(int64_t)(u & 1);
      else if (e.is_signed && bits < 64)
         v = (int64_t)(u ^ (1ULL << (bits - 1))) - (int64_t)(1ULL << (bits - 1));
      if (e.delta && i > 0)
         v = (int64_t)((uint64_t)v + (uint64_t)c->ints[i - 1]);
      c->ints[i] = v;
   }
   release(f, page, 0);
}

/* read_page under a trap: 0 after an error, which leaves the cursor empty. */
static int load(RntCursor *c)
{
   RntFile *f = c->file;
   if (f->err[0] != '\0')
      return 0;
   jmp_buf trap;
   jmp_buf *const outer = f->trap;
   f->trap = &trap;
   if (setjmp(trap) != 0) {
      free_temporaries(f);
      f->trap = outer;
      c->n = c->at = 0;
      return 0;
   }
   read_page(c);
   f->trap = outer;
   return 1;
}

void rnt_seek(RntCursor *c, uint64_t element)
{
   const RntColumn *col = c->col;
   /* The first page that ends after the element. */
   size_t lo = 0;
   size_t hi = col->npages;
   while (lo < hi) {
      const size_t mid = lo + (hi - lo) / 2;
      if (col->pages[mid].first + col->pages[mid].nelem <= element)
         lo = mid + 1;
      else
         hi = mid;
   }
   c->page = lo;
   c->n = c->at = 0;
   c->prev = 0;
   if (lo == col->npages) { /* the end of the column: reading from here fails */
      c->consumed = element;
      c->cluster = UINT32_MAX;
      return;
   }
   const Page *pg = &col->pages[lo];
   if (element < pg->first) {
      fail(c->file, "the pages of a column do not follow each other");
      return;
   }
   c->consumed = pg->first;
   c->cluster = pg->cluster;
   if (!load(c))
      return;
   c->at = (size_t)(element - pg->first);
   if (!encoding(col->type).index)
      return;
   /* Offsets count from the start of the cluster: the one before the element, if in the same cluster. */
   if (c->at > 0)
      c->prev = c->ints[c->at - 1];
   else if (lo > 0 && col->pages[lo - 1].cluster == pg->cluster) {
      rnt_seek(c, element - 1);
      if (c->at < c->n)
         c->prev = c->ints[c->at++];
   }
}

void rnt_seek_items(RntCursor *items, const RntCursor *index)
{
   if (index->cluster == UINT32_MAX) {
      rnt_seek(items, UINT64_MAX);
      return;
   }
   if (index->page == 0 || index->cluster >= items->col->nclusters) {
      fail(items->file, "seek in a collection before its index");
      return;
   }
   rnt_seek(items, items->col->cluster_first[index->cluster] + (uint64_t)index->prev);
}

static int next(RntCursor *c)
{
   while (c->at == c->n) {
      if (!load(c))
         return 0;
   }
   return 1;
}

int64_t rnt_int(RntCursor *c)
{
   if (!next(c))
      return 0;
   if (encoding(c->col->type).size == 1) {
      const uint8_t b = c->bytes[c->at++];
      return encoding(c->col->type).is_signed ? (int64_t)(int8_t)b : (int64_t)b;
   }
   return c->ints[c->at++];
}

uint64_t rnt_len(RntCursor *c)
{
   if (!encoding(c->col->type).index) {
      fail(c->file, "not an index column");
      return 0;
   }
   if (!next(c))
      return 0;
   const int64_t end = c->ints[c->at++];
   if (end < c->prev) {
      fail(c->file, "collection offsets go backwards");
      return 0;
   }
   const uint64_t len = (uint64_t)(end - c->prev);
   c->prev = end;
   return len;
}

void rnt_bytes(RntCursor *c, size_t n, uint8_t *dst)
{
   if (encoding(c->col->type).size != 1)
      fail(c->file, "not a column of bytes");
   while (n > 0) {
      if (c->file->err[0] != '\0' || !next(c)) {
         if (dst)
            memset(dst, 0, n);
         return;
      }
      const size_t k = c->n - c->at < n ? c->n - c->at : n;
      if (dst) {
         memcpy(dst, c->bytes + c->at, k);
         dst += k;
      }
      c->at += k;
      n -= k;
   }
}
