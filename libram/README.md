# libram (proof of concept)

Can a RAM file be read without ROOT? `ramread` writes a RAM file out as SAM,
the same bytes as `ramdump`, with only zstd, the xxHash header and htscodecs,
which htslib already ships. It is C11, for htslib.

- `rntuple.c`, `rntuple.h`: reads the columns of an RNTuple in a ROOT file,
  following the RNTuple binary format specification 1.0. It finds the anchor in
  the ROOT file's top directory, checks the anchor, every envelope and every
  page against its XXH3 checksum, bounds-checks every read, and decodes the
  integer, index, character and byte columns a page at a time.
- `ramread.c`: RAM on top of it: METADATA, the record columns, the tag orders,
  the blocks of names (tok3) and qualities (rANS or fqzcomp), and the PNEXT
  and TLEN predictions. It follows ramtools' RAMNTupleRecord, RAMTags and
  RecordBlocks, so a change to the format needs the same change here.

Build, with htscodecs from pkg-config or from its source:

    make HTSCODECS=<htscodecs source>   # CPPFLAGS=-I<dir> if xxhash.h is elsewhere
    ./ramread -h file.ram > file.sam

`check.sh <ramtools build>/tools ./ramread` writes RAM files with each writer
option and fails if `ramread` prints anything other than what `ramdump` prints.
CI runs it, so a change to the format or to ROOT's encodings shows up there.
Features of RNTuple this reader does not know (feature flags, column types,
compression algorithms) stop it with a message rather than a wrong answer.

Not done yet: region queries, files compressed with anything but ZSTD (the
default), and writing. A damaged file ends the program with a message rather
than returning an error to a caller.
