# Builds without ROOT. Needs htslib, zstd, the xxHash header and htscodecs, from
# pkg-config or from its source: make HTSCODECS=<htscodecs source directory>.
CFLAGS ?= -O2 -g
WARNINGS = -std=c11 -D_POSIX_C_SOURCE=200809L -D_FILE_OFFSET_BITS=64 \
           -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wnull-dereference

ifdef HTSCODECS
# The codecs ramtools' CMake build compiles, without SIMD (an empty config.h).
CODECS_OBJ = fqzcomp_qual.o tokenise_name3.o arith_dynamic.o rANS_static4x16pr.o rANS_static32x16pr.o \
             pack.o rle.o utils.o
CODECS_CFLAGS = -I$(HTSCODECS)
CODECS_LIBS = -lpthread -lm
else
CODECS_CFLAGS = $(shell pkg-config --cflags htscodecs)
CODECS_LIBS = $(shell pkg-config --libs htscodecs)
endif

HTSLIB_CFLAGS = $(shell pkg-config --cflags htslib)
HTSLIB_LIBS = $(shell pkg-config --libs htslib)

ramread: ramread.o ram.o rntuple.o $(CODECS_OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(HTSLIB_LIBS) -lzstd $(CODECS_LIBS)

ramread.o ram.o rntuple.o: %.o: %.c ram.h rntuple.h
	$(CC) $(WARNINGS) $(HTSLIB_CFLAGS) $(CODECS_CFLAGS) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(CODECS_OBJ): %.o: $(HTSCODECS)/htscodecs/%.c
	@mkdir -p codecs && touch codecs/config.h
	$(CC) -Icodecs $(CFLAGS) -c -o $@ $<

clean:
	rm -rf ramread *.o codecs
