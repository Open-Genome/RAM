# Builds without ROOT: needs only zstd and the xxHash header.
CFLAGS ?= -O2 -g
WARNINGS = -std=c11 -D_POSIX_C_SOURCE=200809L -D_FILE_OFFSET_BITS=64 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wnull-dereference

ramread: ramread.c
	$(CC) $(WARNINGS) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(LDFLAGS) -lzstd

clean:
	rm -f ramread
