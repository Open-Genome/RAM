# libram (proof of concept)

Can a RAM file be read without ROOT? `libram` hands back its records as
htslib's `bam1_t`, with only htslib, zstd, the xxHash header and htscodecs,
which htslib already ships. It is C11.

- `rntuple.c`, `rntuple.h`: reads the columns of an RNTuple in a ROOT file,
  following the RNTuple binary format specification 1.0. It finds the anchor in
  the ROOT file's top directory, checks the anchor, every envelope and every
  page against its XXH3 checksum, bounds-checks every read, and decodes the
  integer, index, character and byte columns a page at a time.
- `ram.c`, `ram.h`: RAM on top of it: METADATA, the record columns, the tag
  orders, the blocks of names (tok3) and qualities (rANS or fqzcomp), and the
  PNEXT and TLEN predictions. It follows ramtools' RAMNTupleRecord, RAMTags and
  RecordBlocks, so a change to the format needs the same change here.
- `ramread.c`: writes a RAM file, or a region of it, as SAM or BAM through
  htslib.

The API, as htslib's own readers look:

    RamFile *f = ram_open("file.ram", err, sizeof err);  /* NULL with the reason in err */
    ram_query(f, "chr20:30000000-30100000");              /* optional: a region, as ramdump reads it */
    while ((r = ram_next(f, b)) > 0)                      /* b: a bam1_t; tids index ram_header(f) */
       sam_write1(out, ram_header(f), b);
    if (r < 0) fprintf(stderr, "%s\n", ram_error(f));
    ram_close(f);

Build, with htscodecs from pkg-config or from its source:

    make HTSCODECS=<htscodecs source>   # CPPFLAGS=-I<dir> if xxhash.h is elsewhere
    ./ramread -h file.ram > file.sam
    ./ramread -b file.ram > file.bam
    ./ramread file.ram chr20:30000000-30100000   # a region, as ramdump reads it

The records are the ones htslib makes from the SAM `ramdump` prints, down to
the integer tag types, so `ramread` prints what `samtools view` prints of that
SAM. That is `ramdump`'s output byte for byte unless the input spelled a tag
in a way htslib rewrites, such as `NM:i:007`. A reference the header does not
list cannot be a `bam1_t`, so a record on one is an error.

## samtools

`htslib.sh <empty directory>` builds htslib and samtools 1.22.1 with libram
and `htslib-ram.patch`, so that samtools reads a RAM file as it reads a BAM:

    samtools view -c file.ram chr20:30000000-30100000
    samtools flagstat file.ram
    samtools view -b -o file.bam file.ram

The patch (about 130 lines) makes `hts_open` recognise a ROOT file as RAM and
adds RAM to the switches where htslib picks a reader: header, records, close,
and regions. As for CRAM, the "index" is the open file itself: `ram_query`
needs no index file. Regions follow samtools here, so `chr20:30000000` runs to
the end of chr20. Not done: multi-region iterators (`view -M`), and writing.

`check.sh <ramtools build>/tools ./ramread` writes RAM files with each writer
option and fails if `ramread` prints anything other than what `ramdump` prints,
as SAM, as BAM read back by samtools, and for a region. CI runs it, so a change
to the format or to ROOT's encodings shows up there. Features of RNTuple this
reader does not know (feature flags, column types, compression algorithms)
stop it with a message rather than a wrong answer.

Neither library exits: a damaged file makes `ram_open` return NULL with the
reason, or, met while reading, makes `ram_next` return -1 from then on with
`ram_error` saying why. Underneath, `rnt_open` returns NULL the same way, and
a damaged page makes every cursor return zeros and `rnt_error` say why.

Not done yet: files compressed with anything but ZSTD (the default), and
writing.
