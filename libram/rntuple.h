/* rntuple: read the columns of an RNTuple in a ROOT file without ROOT.

   Follows the RNTuple binary format specification 1.0. Every anchor, envelope
   and page is checked against its XXH3 checksum and every read is bounds
   checked. A damaged file is an error, never a crash or an exit: rnt_open
   returns NULL, and once a cursor meets a damaged page every cursor of the
   file returns zeros and rnt_error says why. An RntFile and its cursors
   belong to one thread. */

#ifndef LIBRAM_RNTUPLE_H
#define LIBRAM_RNTUPLE_H

#include <stddef.h>
#include <stdint.h>

typedef struct RntFile RntFile;
typedef struct RntColumn RntColumn;

/* Reads one column from its first element to its last, a page at a time. */
typedef struct {
   RntFile *file;
   const RntColumn *col;
   size_t page; /* the next page to load */
   uint8_t *bytes;
   int64_t *ints; /* the page's elements, decoded, for columns wider than a byte */
   size_t n, at;  /* elements in the loaded page, and the next one */
   uint32_t cluster;
   int64_t prev;      /* index columns: the offset before the next one, in its cluster */
   uint64_t consumed; /* elements read before the loaded page */
} RntCursor;

/* The RNTuple called name in the ROOT file at path, or NULL with the reason in err. */
RntFile *rnt_open(const char *path, const char *name, char *err, size_t errlen);
void rnt_close(RntFile *f);
uint64_t rnt_entries(const RntFile *f);
/* The first error met reading the file, or NULL. */
const char *rnt_error(const RntFile *f);

/* The which-th column (0 the principal one) of the field at a dotted path such as "record.seq", or NULL
   with an error recorded. */
const RntColumn *rnt_column(RntFile *f, const char *path, unsigned which);

/* A cursor at the first element of col; with a NULL col, one that only returns zeros. */
void rnt_cursor(RntCursor *c, RntFile *f, const RntColumn *col);
void rnt_cursor_free(RntCursor *c);
/* The next value of an integer column. */
int64_t rnt_int(RntCursor *c);
/* The number of items in the next collection, from its index column. */
uint64_t rnt_len(RntCursor *c);
/* The next n elements of a column of bytes or characters, into dst (NULL skips them). */
void rnt_bytes(RntCursor *c, size_t n, uint8_t *dst);
/* Moves the cursor to the element at that index in its column: an entry, for a top-level field. */
void rnt_seek(RntCursor *c, uint64_t element);
/* Moves items to the first item of the collection that index, just moved with rnt_seek, reads next. */
void rnt_seek_items(RntCursor *items, const RntCursor *index);

#endif
