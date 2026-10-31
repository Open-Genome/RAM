/* ram: read a RAM file without ROOT, a record at a time, as htslib's bam1_t.

   A damaged file is an error, never a crash or an exit: ram_open returns NULL,
   ram_next returns -1 from then on, and ram_error says why. A RamFile belongs
   to one thread. */

#ifndef LIBRAM_RAM_H
#define LIBRAM_RAM_H

#include <htslib/sam.h>

#include <stddef.h>

typedef struct RamFile RamFile;

/* The RAM file at path, set to read every record, or NULL with the reason in err. */
RamFile *ram_open(const char *path, char *err, size_t errlen);
void ram_close(RamFile *f);
/* The file's header, which the file owns. Record tids index it. */
sam_hdr_t *ram_header(RamFile *f);
/* From the next ram_next on, the records ramdump prints for region: rname[:start[-end]], 1-based and
   inclusive, where rname:pos is that base alone. NULL or "*" means every record. Returns 0, or -1 for a
   region the file cannot have, with ram_error saying why; the file stays usable. */
int ram_query(RamFile *f, const char *region);
/* Reads the next record into b: 1, or 0 at the end, or -1 on an error. */
int ram_next(RamFile *f, bam1_t *b);
/* Why the last call failed, or NULL. */
const char *ram_error(const RamFile *f);

#endif
