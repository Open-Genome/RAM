# libram (proof of concept)

Can a RAM file be read without ROOT? `ramread` answers for the container:
it is about 550 lines of C11 that need only zstd and the xxHash header, and
print one integer field of a RAM file.

It follows the RNTuple binary format specification 1.0: it finds the
`ROOT::RNTuple` anchor of `RAM` in the ROOT file's top directory, reads the
header, footer and page list envelopes, and decodes the pages of the field.
Every envelope, page and the anchor is checked against its XXH3 checksum, and
every read is bounds-checked, since the file is untrusted.

    make                       # CPPFLAGS=-I<dir> if xxhash.h is elsewhere
    ./ramread file.ram record.pos

It reads the integer fields of the record (`pos`, `flag`, `mapq`, `refid`,
`tag_nm`, ...) of files written with ZSTD, the default. Strings, vectors,
the name and quality blocks, METADATA and writing are not done yet.
