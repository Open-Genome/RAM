/* ramread: print one integer column of a RAM file without ROOT.

   A proof that the RNTuple binary format (specification 1.0) can be read with
   zstd and xxHash alone: it finds the anchor in the ROOT file, checks every
   envelope and page against its checksum and decodes the pages of one field.

   Usage: ramread <file.ram> <field>      e.g. ramread x.ram record.pos */

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

_Noreturn static void die(const char *fmt, ...)
{
   va_list ap;
   va_start(ap, fmt);
   fputs("ramread: ", stderr);
   vfprintf(stderr, fmt, ap);
   fputc('\n', stderr);
   va_end(ap);
   exit(1);
}

static void *xmalloc(size_t n)
{
   void *p = malloc(n ? n : 1);
   if (p == NULL)
      die("out of memory");
   return p;
}

/* A byte range being parsed. Every read is checked against its end: the file is untrusted. */
typedef struct {
   const uint8_t *p;
   size_t n, off;
} Cursor;

static const uint8_t *take(Cursor *c, size_t n)
{
   if (n > c->n - c->off)
      die("truncated record");
   const uint8_t *q = c->p + c->off;
   c->off += n;
   return q;
}

static void seek(Cursor *c, uint64_t off)
{
   if (off > c->n)
      die("offset past the end of a record");
   c->off = (size_t)off;
}

static uint64_t le(Cursor *c, size_t n)
{
   const uint8_t *q = take(c, n);
   uint64_t v = 0;
   for (size_t i = n; i-- > 0;)
      v = (v << 8) | q[i];
   return v;
}

static uint64_t be(Cursor *c, size_t n)
{
   const uint8_t *q = take(c, n);
   uint64_t v = 0;
   for (size_t i = 0; i < n; ++i)
      v = (v << 8) | q[i];
   return v;
}

static uint64_t le64_at(const uint8_t *q)
{
   Cursor c = {q, 8, 0};
   return le(&c, 8);
}

/* A ROOT TString: one length byte, or 255 followed by a 32-bit big-endian length. */
static void tstring(Cursor *c, const uint8_t **s, size_t *len)
{
   size_t n = (size_t)be(c, 1);
   if (n == 255)
      n = (size_t)be(c, 4);
   *s = take(c, n);
   *len = n;
}

/* An RNTuple string: a 32-bit little-endian length and the bytes. */
static char *ntstring(Cursor *c)
{
   const size_t n = (size_t)le(c, 4);
   const uint8_t *s = take(c, n);
   char *out = xmalloc(n + 1);
   memcpy(out, s, n);
   out[n] = '\0';
   return out;
}

/* Opens a frame and returns the offset where it ends; a list frame also gives its number of items. */
static size_t frame(Cursor *c, int want_list, uint32_t *nitems)
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
   FILE *f;
   uint64_t size;
} File;

static uint8_t *read_at(const File *f, uint64_t off, uint64_t n)
{
   if (off > f->size || n > f->size - off)
      die("range %" PRIu64 "+%" PRIu64 " is past the end of the file", off, n);
   uint8_t *buf = xmalloc((size_t)n);
   if (fseeko(f->f, (off_t)off, SEEK_SET) != 0 || fread(buf, 1, (size_t)n, f->f) != n)
      die("cannot read %" PRIu64 " bytes at %" PRIu64, n, off);
   return buf;
}

/* Undoes ROOT's compression blocks: chunks of a 9-byte header and the compressed bytes. */
static uint8_t *unzip(uint8_t *in, size_t nin, size_t nout)
{
   if (nin == nout)
      return in;
   uint8_t *out = xmalloc(nout);
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
envelope(const File *f, uint64_t at, uint64_t nbytes, uint64_t len, unsigned type, Cursor *c, uint64_t *checksum)
{
   if (len < 16 || len > SIZE_MAX)
      die("bad envelope length");
   uint8_t *buf = unzip(read_at(f, at, nbytes), (size_t)nbytes, (size_t)len);
   const uint64_t ck = le64_at(buf + len - 8);
   if (XXH3_64bits(buf, (size_t)len - 8) != ck)
      die("envelope %u checksum mismatch", type);
   *c = (Cursor){buf, (size_t)len - 8, 0};
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

static Locator locator(Cursor *c)
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

static void check_features(Cursor *c)
{
   if (le(c, 8) != 0)
      die("the file uses RNTuple features this reader does not know");
}

enum {
   MAX_FIELDS = 4096
};

typedef struct {
   uint32_t parent;
   char *name;
} Field;

typedef struct {
   uint16_t type;
   uint16_t repr;
   uint32_t field;
} Column;

typedef struct {
   Field fields[MAX_FIELDS];
   Column columns[MAX_FIELDS];
   uint32_t nfields, ncolumns;
} Schema;

/* The four list frames of a schema description, in the header or its extension in the footer. */
static void schema(Cursor *c, Schema *s)
{
   uint32_t n = 0;
   size_t end = frame(c, 1, &n);
   for (uint32_t i = 0; i < n; ++i) {
      if (s->nfields == MAX_FIELDS)
         die("too many fields");
      const size_t fend = frame(c, 0, NULL);
      take(c, 8); /* field and type version */
      Field *fd = &s->fields[s->nfields++];
      fd->parent = (uint32_t)le(c, 4);
      take(c, 4); /* structural role, flags */
      fd->name = ntstring(c);
      c->off = fend;
   }
   c->off = end;

   end = frame(c, 1, &n);
   for (uint32_t i = 0; i < n; ++i) {
      if (s->ncolumns == MAX_FIELDS)
         die("too many columns");
      const size_t cend = frame(c, 0, NULL);
      Column *col = &s->columns[s->ncolumns++];
      col->type = (uint16_t)le(c, 2);
      take(c, 2); /* bits on storage */
      col->field = (uint32_t)le(c, 4);
      take(c, 2); /* flags */
      col->repr = (uint16_t)le(c, 2);
      c->off = cend;
   }
   c->off = end;

   for (int k = 0; k < 2; ++k) { /* alias columns, extra type information */
      end = frame(c, 1, &n);
      c->off = end;
   }
}

/* The field ID of a dotted path such as "record.pos", or -1. */
static int64_t find_field(const Schema *s, const char *path)
{
   int64_t parent = -1;
   while (*path) {
      const char *dot = strchr(path, '.');
      const size_t len = dot ? (size_t)(dot - path) : strlen(path);
      int64_t found = -1;
      for (uint32_t i = 0; i < s->nfields && found < 0; ++i) {
         const int top = s->fields[i].parent == i;
         const int under = parent < 0 ? top : (!top && s->fields[i].parent == (uint64_t)parent);
         if (under && strlen(s->fields[i].name) == len && memcmp(s->fields[i].name, path, len) == 0)
            found = i;
      }
      if (found < 0)
         return -1;
      parent = found;
      path += len + (dot != NULL);
   }
   return parent;
}

/* Element size, and whether the column is split and zigzag encoded, for the integer column types. */
static int int_column(uint16_t type, size_t *size, int *split, int *is_signed)
{
   switch (type) {
   case 0x03: *size = 1, *split = 0, *is_signed = 1; return 1;
   case 0x04: *size = 1, *split = 0, *is_signed = 0; return 1;
   case 0x05: *size = 2, *split = 0, *is_signed = 1; return 1;
   case 0x06: *size = 2, *split = 0, *is_signed = 0; return 1;
   case 0x07: *size = 4, *split = 0, *is_signed = 1; return 1;
   case 0x08: *size = 4, *split = 0, *is_signed = 0; return 1;
   case 0x09: *size = 8, *split = 0, *is_signed = 1; return 1;
   case 0x0A: *size = 8, *split = 0, *is_signed = 0; return 1;
   case 0x11: *size = 2, *split = 1, *is_signed = 1; return 1;
   case 0x12: *size = 2, *split = 1, *is_signed = 0; return 1;
   case 0x13: *size = 4, *split = 1, *is_signed = 1; return 1;
   case 0x14: *size = 4, *split = 1, *is_signed = 0; return 1;
   case 0x15: *size = 8, *split = 1, *is_signed = 1; return 1;
   case 0x16: *size = 8, *split = 1, *is_signed = 0; return 1;
   default: return 0;
   }
}

static void print_page(const File *f, int32_t nelem_ck, Locator loc, const Column *col)
{
   size_t size = 0;
   int split = 0;
   int is_signed = 0;
   if (!int_column(col->type, &size, &split, &is_signed))
      die("not an integer column");
   const int has_ck = nelem_ck < 0;
   const uint64_t nelem = has_ck ? (uint64_t)0 - (uint64_t)(int64_t)nelem_ck : (uint64_t)nelem_ck;
   if (nelem > (1U << 30) / size || loc.size > UINT64_MAX - 8)
      die("bad page size");
   const size_t nout = (size_t)nelem * size;

   uint8_t *raw = read_at(f, loc.offset, loc.size + (has_ck ? 8 : 0));
   if (has_ck && XXH3_64bits(raw, (size_t)loc.size) != le64_at(raw + loc.size))
      die("page checksum mismatch at %" PRIu64, loc.offset);
   uint8_t *page = unzip(raw, (size_t)loc.size, nout);

   for (size_t i = 0; i < nelem; ++i) {
      uint64_t u = 0;
      for (size_t b = size; b-- > 0;)
         u = (u << 8) | page[split ? b * (size_t)nelem + i : i * size + b];
      if (!is_signed) {
         printf("%" PRIu64 "\n", u);
         continue;
      }
      int64_t v = 0;
      if (split)
         v = (int64_t)(u >> 1) ^ -(int64_t)(u & 1);
      else
         v = size == 8 ? (int64_t)u : (int64_t)(u ^ (1ULL << (8 * size - 1))) - (int64_t)(1ULL << (8 * size - 1));
      printf("%" PRId64 "\n", v);
   }
   free(page);
}

/* Finds the key of the RNTuple called name in the top directory and returns the offset of its object. */
static uint64_t find_anchor(const File *f, const char *name, uint64_t *objlen)
{
   uint8_t *hb = read_at(f, 0, f->size < 64 ? f->size : 64);
   Cursor c = {hb, f->size < 64 ? (size_t)f->size : 64, 0};
   if (memcmp(take(&c, 4), "root", 4) != 0)
      die("not a ROOT file");
   take(&c, 4); /* version */
   const uint64_t begin = be(&c, 4);
   free(hb);
   if (begin >= f->size)
      die("bad file header");

   /* The top directory: its key, the file's name and title, then the streamed TDirectory. */
   uint8_t *db = read_at(f, begin, f->size - begin < 512 ? f->size - begin : 512);
   c = (Cursor){db, f->size - begin < 512 ? (size_t)(f->size - begin) : 512, 0};
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
   c = (Cursor){kb, (size_t)nbytes_keys, 0};
   take(&c, 14);
   seek(&c, be(&c, 2));
   const uint64_t nkeys = be(&c, 4);
   uint64_t seek = 0;
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
         seek = seek_key;
         *objlen = len_obj;
      }
   }
   free(kb);
   if (seek == 0)
      die("no RNTuple called %s", name);

   /* The object follows the key's own header. */
   uint8_t *ab = read_at(f, seek, 18);
   c = (Cursor){ab, 18, 0};
   const uint64_t nbytes = be(&c, 4);
   take(&c, 10);
   const uint64_t keylen = be(&c, 2);
   free(ab);
   if (nbytes < keylen || nbytes - keylen != *objlen)
      die("compressed anchors are not supported");
   return seek + keylen;
}

int main(int argc, char **argv)
{
   if (argc != 3) {
      fprintf(stderr, "Usage: ramread <file.ram> <field>\n");
      return 2;
   }
   File f = {fopen(argv[1], "rb"), 0};
   struct stat st;
   if (f.f == NULL || fstat(fileno(f.f), &st) != 0)
      die("cannot open %s", argv[1]);
   f.size = (uint64_t)st.st_size;

   /* The anchor: big-endian, with an XXH3 over everything after byte count and class version. */
   uint64_t objlen = 0;
   const uint64_t at = find_anchor(&f, "RAM", &objlen);
   if (objlen < 78)
      die("anchor too short");
   uint8_t *ab = read_at(&f, at, objlen);
   Cursor c = {ab, (size_t)objlen, 0};
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

   static Schema s;
   uint64_t header_ck = 0;
   uint8_t *hb = envelope(&f, seek_header, nbytes_header, len_header, 1, &c, &header_ck);
   check_features(&c);
   for (int i = 0; i < 3; ++i)
      free(ntstring(&c)); /* name, description, writer */
   schema(&c, &s);
   free(hb);

   uint8_t *fb = envelope(&f, seek_footer, nbytes_footer, len_footer, 2, &c, NULL);
   check_features(&c);
   if (le(&c, 8) != header_ck)
      die("footer belongs to another header");
   const size_t ext_end = frame(&c, 0, NULL);
   if (c.off < ext_end)
      schema(&c, &s);
   c.off = ext_end;

   const int64_t field = find_field(&s, argv[2]);
   if (field < 0)
      die("no field %s", argv[2]);
   int64_t column = -1;
   for (uint32_t i = 0; i < s.ncolumns; ++i) {
      if (s.columns[i].field == (uint64_t)field) {
         if (column >= 0 || s.columns[i].repr != 0)
            die("%s is not a single integer column", argv[2]);
         column = i;
      }
   }
   size_t size = 0;
   int split = 0;
   int is_signed = 0;
   if (column < 0 || !int_column(s.columns[column].type, &size, &split, &is_signed))
      die("%s is not a single integer column", argv[2]);

   uint32_t ngroups = 0;
   const size_t groups_end = frame(&c, 1, &ngroups);
   for (uint32_t g = 0; g < ngroups; ++g) {
      const size_t gend = frame(&c, 0, NULL);
      take(&c, 20); /* min entry, entry span, number of clusters */
      const uint64_t len_pl = le(&c, 8);
      const Locator pl = locator(&c);
      c.off = gend;

      Cursor p;
      uint8_t *pb = envelope(&f, pl.offset, pl.size, len_pl, 3, &p, NULL);
      if (le(&p, 8) != header_ck)
         die("page list belongs to another header");
      uint32_t n = 0;
      size_t end = frame(&p, 1, &n);
      for (uint32_t i = 0; i < n; ++i) {
         const size_t send = frame(&p, 0, NULL);
         take(&p, 8);
         if (le(&p, 8) >> 56 & 1)
            die("sharded clusters are not supported");
         p.off = send;
      }
      p.off = end;

      uint32_t nclusters = 0;
      const size_t clusters_end = frame(&p, 1, &nclusters);
      for (uint32_t k = 0; k < nclusters; ++k) {
         uint32_t ncols = 0;
         const size_t cols_end = frame(&p, 1, &ncols);
         for (uint32_t col = 0; col < ncols; ++col) {
            uint32_t npages = 0;
            const size_t pages_end = frame(&p, 1, &npages);
            for (uint32_t pg = 0; pg < npages; ++pg) {
               const int32_t nelem = (int32_t)le(&p, 4);
               const Locator loc = locator(&p);
               if (col == column)
                  print_page(&f, nelem, loc, &s.columns[column]);
            }
            p.off = pages_end;
         }
         p.off = cols_end;
      }
      p.off = clusters_end;
      free(pb);
   }
   c.off = groups_end;
   free(fb);
   for (uint32_t i = 0; i < s.nfields; ++i)
      free(s.fields[i].name);
   fclose(f.f);
   return 0;
}
