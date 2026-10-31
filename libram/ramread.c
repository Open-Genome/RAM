/* ramread: write a RAM file, or a region of it, out as SAM, as ramdump does, or as BAM, without ROOT.

   Usage: ramread [-h] [-b] <file.ram> [region] */

#include "ram.h"

#include <htslib/sam.h>

#include <stdio.h>
#include <string.h>

int main(int argc, char **argv)
{
   int header = 0;
   int bam = 0;
   int usage = argc < 2;
   const char *path = NULL;
   const char *region = NULL;
   for (int i = 1; i < argc; ++i) {
      if (strcmp(argv[i], "-h") == 0)
         header = 1;
      else if (strcmp(argv[i], "-b") == 0)
         bam = 1;
      else if (path == NULL && argv[i][0] != '-')
         path = argv[i];
      else if (region == NULL && argv[i][0] != '-')
         region = argv[i];
      else
         usage = 1;
   }
   if (usage || path == NULL) {
      fputs("Usage: ramread [-h] [-b] <file.ram> [region]\n"
            "-h writes the header and -b writes BAM, which always has it.\n"
            "Region is rname[:start[-end]], 1-based and inclusive, as in samtools.\n",
            stderr);
      return 2;
   }

   char err[256];
   RamFile *f = ram_open(path, err, sizeof err);
   if (f == NULL) {
      fprintf(stderr, "ramread: %s\n", err);
      return 1;
   }
   int status = 1;
   const char *why = "cannot write the output";
   samFile *out = NULL;
   bam1_t *b = bam_init1();
   if (region != NULL && ram_query(f, region) < 0)
      why = ram_error(f);
   else
      out = sam_open("-", bam ? "wb" : "w");
   if (out != NULL && b != NULL && (!(header || bam) || sam_hdr_write(out, ram_header(f)) == 0)) {
      int got = 0;
      while ((got = ram_next(f, b)) > 0 && sam_write1(out, ram_header(f), b) >= 0) {
      }
      if (got < 0)
         why = ram_error(f);
      status = got != 0;
   }
   bam_destroy1(b);
   if (out != NULL && sam_close(out) < 0)
      status = 1;
   if (status)
      fprintf(stderr, "ramread: %s\n", why);
   ram_close(f);
   return status;
}
