#include "rntuple.h"

#define XXH_INLINE_ALL
#include <xxhash.h>
#include <zstd.h>

#include <inttypes.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

void rnt_die(const char *fmt, ...)
{
   va_list ap;
   va_start(ap, fmt);
   fputs("ramread: ", stderr);
   vfprintf(stderr, fmt, ap);
   fputc('\n', stderr);
   va_end(ap);
   exit(1);
}

#define die rnt_die

static void *xrealloc(void *p, size_t n)
{
   void *q = realloc(p, n ? n : 1);
   if (q == NULL)
      die("out of memory");
   return q;
}

/* A byte range being parsed. Every read is checked against its end: the file is untrusted. */
typedef struct {
   const uint8_t *p;
   size_t n, off;
} Buf;

static const uint8_t *take(Buf *c, size_t n)
{
   if (n > c->n - c->off)
      die("truncated record");
   const uint8_t *q = c->p + c->off;
   c->off += n;
   return q;
}

static void seek(Buf *c, uint64_t off)
{
   if (off > c->n)
      die("offset past the end of a record");
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
   Buf c = {q, 8, 0};
   return le(&c, 8);
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
   char *out = xrealloc(NULL, n + 1);
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
      die("expected a %s frame", want_list ? "list" : "record");
   if (is_list)
      *nitems = (uint32_t)le(c, 4);
   if (len > c->n - start || start + len < c->off)
      die("bad frame size");
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
};

struct RntFile {
   FILE *f;
   uint64_t size;
   char **names;
   uint32_t *parents;
   uint32_t nfields;
   RntColumn *columns;
   uint32_t ncolumns;
   uint64_t entries;
};

static uint8_t *read_at(const RntFile *f, uint64_t off, uint64_t n)
{
   if (off > f->size || n > f->size - off)
      die("range %" PRIu64 "+%" PRIu64 " is past the end of the file", off, n);
   uint8_t *buf = xrealloc(NULL, (size_t)n);
   if (fseeko(f->f, (off_t)off, SEEK_SET) != 0 || fread(buf, 1, (size_t)n, f->f) != n)
      die("cannot read %" PRIu64 " bytes at %" PRIu64, n, off);
   return buf;
}

/* Undoes ROOT's compression blocks: chunks of a 9-byte header and the compressed bytes. */
static uint8_t *unzip(uint8_t *in, size_t nin, size_t nout)
{
   if (nin == nout)
      return in;
   uint8_t *out = xrealloc(NULL, nout);
   size_t i = 0;
   size_t o = 0;
   while (o < nout) {
      if (nin - i < 9)
         die("truncated compression block");
      const uint8_t *h = in + i;
      const size_t csize = (size_t)h[3] | (size_t)h[4] << 8 | (size_t)h[5] << 16;
      const size_t usize = (size_t)h[6] | (size_t)h[7] << 8 | (size_t)h[8] << 16;
      if (h[0] != 'Z' || h[1] != 'S')
         die("compression %c%c is not supported, only zstd", h[0], h[1]);
      if (csize > nin - i - 9 || usize > nout - o)
         die("bad compression block");
      if (ZSTD_decompress(out + o, usize, h + 9, csize) != usize)
         die("zstd failed");
      i += 9 + csize;
      o += usize;
   }
   if (i != nin)
      die("bad compression block");
   free(in);
   return out;
}

/* Reads an envelope, checks its type and checksum, and leaves the cursor on its payload. */
static uint8_t *
envelope(const RntFile *f, uint64_t at, uint64_t nbytes, uint64_t len, unsigned type, Buf *c, uint64_t *checksum)
{
   if (len < 16 || len > SIZE_MAX)
      die("bad envelope length");
   uint8_t *buf = unzip(read_at(f, at, nbytes), (size_t)nbytes, (size_t)len);
   const uint64_t ck = le64_at(buf + len - 8);
   if (XXH3_64bits(buf, (size_t)len - 8) != ck)
      die("envelope %u checksum mismatch", type);
   *c = (Buf){buf, (size_t)len - 8, 0};
   const uint64_t id = le(c, 8);
   if ((id & 0xffff) != type || id >> 16 != len)
      die("expected envelope %u of %" PRIu64 " bytes", type, len);
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
      die("locator type %d is not supported", type);
   const uint64_t size = le(c, 8);
   return (Locator){size, le(c, 8)};
}

static void check_features(Buf *c)
{
   if (le(c, 8) != 0)
      die("the file uses RNTuple features this reader does not know");
}

/* The four list frames of a schema description, in the header or its extension in the footer. */
static void schema(Buf *c, RntFile *f)
{
   uint32_t n = 0;
   size_t end = frame(c, 1, &n);
   if (n > c->n)
      die("bad field count");
   f->names = xrealloc(f->names, (f->nfields + n) * sizeof *f->names);
   f->parents = xrealloc(f->parents, (f->nfields + n) * sizeof *f->parents);
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
      die("bad column count");
   f->columns = xrealloc(f->columns, (f->ncolumns + n) * sizeof *f->columns);
   for (uint32_t i = 0; i < n; ++i) {
      const size_t cend = frame(c, 0, NULL);
      RntColumn *col = &f->columns[f->ncolumns++];
      memset(col, 0, sizeof *col);
      col->type = (uint16_t)le(c, 2);
      take(c, 2); /* bits on storage */
      col->field = (uint32_t)le(c, 4);
      if (le(c, 2) & 1)
         die("deferred columns are not supported");
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
static uint64_t find_anchor(const RntFile *f, const char *name, uint64_t *objlen)
{
   uint8_t *hb = read_at(f, 0, f->size < 64 ? f->size : 64);
   Buf c = {hb, f->size < 64 ? (size_t)f->size : 64, 0};
   if (memcmp(take(&c, 4), "root", 4) != 0)
      die("not a ROOT file");
   take(&c, 4); /* version */
   const uint64_t begin = be(&c, 4);
   free(hb);
   if (begin >= f->size)
      die("bad file header");

   /* The top directory: its key, the file's name and title, then the streamed TDirectory. */
   uint8_t *db = read_at(f, begin, f->size - begin < 512 ? f->size - begin : 512);
   c = (Buf){db, f->size - begin < 512 ? (size_t)(f->size - begin) : 512, 0};
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
   free(db);

   /* The key list: a key of its own, the number of keys, then each key's header. */
   uint8_t *kb = read_at(f, seek_keys, nbytes_keys);
   c = (Buf){kb, (size_t)nbytes_keys, 0};
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
   free(kb);
   if (at == 0)
      die("no RNTuple called %s", name);

   /* The object follows the key's own header. */
   uint8_t *ab = read_at(f, at, 18);
   c = (Buf){ab, 18, 0};
   const uint64_t nbytes = be(&c, 4);
   take(&c, 10);
   const uint64_t keylen = be(&c, 2);
   free(ab);
   if (nbytes < keylen || nbytes - keylen != *objlen)
      die("compressed anchors are not supported");
   return at + keylen;
}

/* Adds the pages of every column of the clusters in one page list envelope. */
static void page_list(RntFile *f, Buf *p, uint64_t header_ck, uint32_t first_cluster)
{
   if (le(p, 8) != header_ck)
      die("page list belongs to another header");
   uint32_t n = 0;
   size_t end = frame(p, 1, &n);
   for (uint32_t i = 0; i < n; ++i) {
      const size_t send = frame(p, 0, NULL);
      take(p, 8); /* first entry */
      const uint64_t entries = le(p, 8);
      if (entries >> 56 & 1)
         die("sharded clusters are not supported");
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
         die("a cluster has %u columns, the schema %u", ncols, f->ncolumns);
      for (uint32_t i = 0; i < ncols; ++i) {
         RntColumn *col = &f->columns[i];
         uint32_t npages = 0;
         const size_t pages_end = frame(p, 1, &npages);
         if (npages > p->n)
            die("bad page count");
         col->pages = xrealloc(col->pages, (col->npages + npages) * sizeof *col->pages);
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
            die("suppressed columns are not supported");
         for (uint32_t pg = 0; pg < npages; ++pg)
            first[pg].first += (uint64_t)offset;
         p->off = pages_end;
      }
      p->off = cols_end;
   }
   p->off = clusters_end;
}

RntFile *rnt_open(const char *path, const char *name)
{
   RntFile *f = xrealloc(NULL, sizeof *f);
   memset(f, 0, sizeof *f);
   f->f = fopen(path, "rb");
   struct stat st;
   if (f->f == NULL || fstat(fileno(f->f), &st) != 0)
      die("cannot open %s", path);
   f->size = (uint64_t)st.st_size;

   /* The anchor: big-endian, with an XXH3 over everything after byte count and class version. */
   uint64_t objlen = 0;
   const uint64_t at = find_anchor(f, name, &objlen);
   if (objlen < 78)
      die("anchor too short");
   uint8_t *ab = read_at(f, at, objlen);
   Buf c = {ab, (size_t)objlen, 0};
   seek(&c, objlen - 8);
   if (XXH3_64bits(ab + 6, (size_t)objlen - 14) != be(&c, 8))
      die("anchor checksum mismatch");
   seek(&c, 6);
   if (be(&c, 2) != 1)
      die("RNTuple epoch is not 1");
   take(&c, 6); /* major, minor, patch */
   const uint64_t seek_header = be(&c, 8);
   const uint64_t nbytes_header = be(&c, 8);
   const uint64_t len_header = be(&c, 8);
   const uint64_t seek_footer = be(&c, 8);
   const uint64_t nbytes_footer = be(&c, 8);
   const uint64_t len_footer = be(&c, 8);
   free(ab);

   uint64_t header_ck = 0;
   uint8_t *hb = envelope(f, seek_header, nbytes_header, len_header, 1, &c, &header_ck);
   check_features(&c);
   for (int i = 0; i < 3; ++i)
      free(ntstring(&c)); /* name, description, writer */
   schema(&c, f);
   free(hb);

   uint8_t *fb = envelope(f, seek_footer, nbytes_footer, len_footer, 2, &c, NULL);
   check_features(&c);
   if (le(&c, 8) != header_ck)
      die("footer belongs to another header");
   const size_t ext_end = frame(&c, 0, NULL);
   if (c.off < ext_end)
      schema(&c, f);
   c.off = ext_end;
   for (uint32_t i = 0; i < f->ncolumns; ++i) {
      if (f->columns[i].field >= f->nfields || f->columns[i].repr != 0)
         die("column %u is not supported", i);
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
      free(pb);
      clusters += nclusters;
   }
   c.off = groups_end;
   free(fb);
   return f;
}

void rnt_close(RntFile *f)
{
   for (uint32_t i = 0; i < f->nfields; ++i)
      free(f->names[i]);
   for (uint32_t i = 0; i < f->ncolumns; ++i)
      free(f->columns[i].pages);
   free(f->names);
   free(f->parents);
   free(f->columns);
   fclose(f->f);
   free(f);
}

uint64_t rnt_entries(const RntFile *f)
{
   return f->entries;
}

const RntColumn *rnt_column(const RntFile *f, const char *path, unsigned which)
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
      if (found < 0)
         die("no field %s", path);
      parent = found;
      p += len + (dot != NULL);
   }
   for (uint32_t i = 0; i < f->ncolumns; ++i) {
      if (f->columns[i].field == (uint64_t)parent && which-- == 0)
         return &f->columns[i];
   }
   die("field %s has no such column", path);
}

/* How the column type stores an element: its size, and whether it is split, zigzag, delta, signed or an index. */
typedef struct {
   size_t size;
   int split, zigzag, delta, is_signed, index;
} Encoding;

static Encoding encoding(uint16_t type)
{
   switch (type) {
   case 0x01:
   case 0x02:
   case 0x04: return (Encoding){1, 0, 0, 0, 0, 0};
   case 0x03: return (Encoding){1, 0, 0, 0, 1, 0};
   case 0x05: return (Encoding){2, 0, 0, 0, 1, 0};
   case 0x06: return (Encoding){2, 0, 0, 0, 0, 0};
   case 0x07: return (Encoding){4, 0, 0, 0, 1, 0};
   case 0x08: return (Encoding){4, 0, 0, 0, 0, 0};
   case 0x09: return (Encoding){8, 0, 0, 0, 1, 0};
   case 0x0A: return (Encoding){8, 0, 0, 0, 0, 0};
   case 0x0E: return (Encoding){4, 0, 0, 0, 0, 1};
   case 0x0F: return (Encoding){8, 0, 0, 0, 0, 1};
   case 0x11: return (Encoding){2, 1, 1, 0, 1, 0};
   case 0x12: return (Encoding){2, 1, 0, 0, 0, 0};
   case 0x13: return (Encoding){4, 1, 1, 0, 1, 0};
   case 0x14: return (Encoding){4, 1, 0, 0, 0, 0};
   case 0x15: return (Encoding){8, 1, 1, 0, 1, 0};
   case 0x16: return (Encoding){8, 1, 0, 0, 0, 0};
   case 0x1A: return (Encoding){4, 1, 0, 1, 0, 1};
   case 0x1B: return (Encoding){8, 1, 0, 1, 0, 1};
   default: die("column type 0x%02x is not supported", type);
   }
}

void rnt_cursor(RntCursor *c, const RntFile *f, const RntColumn *col)
{
   memset(c, 0, sizeof *c);
   c->file = f;
   c->col = col;
   encoding(col->type);
}

void rnt_cursor_free(RntCursor *c)
{
   free(c->bytes);
   free(c->ints);
   memset(c, 0, sizeof *c);
}

/* Loads the next page of the cursor's column and decodes its elements. */
static void load(RntCursor *c)
{
   const RntColumn *col = c->col;
   if (c->page == col->npages)
      die("a column ends before its last entry");
   const Page *pg = &col->pages[c->page];
   const Encoding e = encoding(col->type);
   c->consumed += c->n;
   if (pg->first != c->consumed)
      die("the pages of a column do not follow each other");
   if (c->page == 0 || pg->cluster != c->cluster)
      c->prev = 0;
   c->cluster = pg->cluster;
   c->page++;

   if (pg->nelem > (1U << 30) / e.size || pg->size > UINT64_MAX - 8)
      die("bad page size");
   const size_t n = pg->nelem;
   uint8_t *raw = read_at(c->file, pg->offset, pg->size + (pg->checksum ? 8 : 0));
   if (pg->checksum && XXH3_64bits(raw, (size_t)pg->size) != le64_at(raw + pg->size))
      die("page checksum mismatch at %" PRIu64, pg->offset);
   uint8_t *page = unzip(raw, (size_t)pg->size, n * e.size);
   c->n = n;
   c->at = 0;
   if (e.size == 1) {
      free(c->bytes);
      c->bytes = page;
      return;
   }
   c->ints = xrealloc(c->ints, (n ? n : 1) * sizeof *c->ints);
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
   free(page);
}

static void next(RntCursor *c)
{
   while (c->at == c->n)
      load(c);
}

int64_t rnt_int(RntCursor *c)
{
   next(c);
   if (encoding(c->col->type).size == 1) {
      const uint8_t b = c->bytes[c->at++];
      return encoding(c->col->type).is_signed ? (int64_t)(int8_t)b : (int64_t)b;
   }
   return c->ints[c->at++];
}

uint64_t rnt_len(RntCursor *c)
{
   if (!encoding(c->col->type).index)
      die("not an index column");
   next(c);
   const int64_t end = c->ints[c->at++];
   if (end < c->prev)
      die("collection offsets go backwards");
   const uint64_t len = (uint64_t)(end - c->prev);
   c->prev = end;
   return len;
}

void rnt_bytes(RntCursor *c, size_t n, uint8_t *dst)
{
   if (encoding(c->col->type).size != 1)
      die("not a column of bytes");
   while (n > 0) {
      next(c);
      const size_t k = c->n - c->at < n ? c->n - c->at : n;
      if (dst) {
         memcpy(dst, c->bytes + c->at, k);
         dst += k;
      }
      c->at += k;
      n -= k;
   }
}
