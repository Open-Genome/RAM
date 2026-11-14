# RAM

RAM (ROOT Alignment/Map) is a file format for aligned sequencing reads: the
same data as a SAM or BAM file, stored column by column in RNTuple, the
columnar format CERN's ROOT uses for petabytes of physics data. A RAM file

- gives back every record byte for byte, header and tags included;
- needs no reference sequence to be read;
- answers region queries such as `chr13:32889611-32973805`;
- is less than half the size of the BAM of the same reads, and close to a
  CRAM 3.1 that also reads without a reference.

On 5 Mb of chr20 (see [Benchmarks](#benchmarks)):

| | BAM | CRAM 3.1, `embed_ref=2` | RAM |
|---|---:|---:|---:|
| Illumina short reads (HG00154, 30x) | 152.4 MB | 69.1 MB | 70.9 MB |
| PacBio Revio HiFi reads (HG002) | 101.2 MB | 41.9 MB | 44.8 MB |

This repository holds the format and its tools:

| Tool | What it does |
|------|--------------|
| `samtoramntuple` | converts a SAM file to RAM |
| `bamtoramntuple` | converts a BAM file to RAM |
| `ramdump` | writes a RAM file, or a region of it, out as SAM; takes the `samtools view` options |
| `ramntupleview` | counts the records in a region and reports the query time |
| `ramstats` | `idxstats` and `flagstat`, printed as `samtools idxstats` and `samtools flagstat` print them, and `size`, the bytes each field takes |

Every result can be checked against samtools: `ramdump` writes SAM, and
`ramstats` prints what samtools prints.

## How RAM stores reads

A RAM file is a ROOT file holding one RNTuple, `RAM`, and its METADATA.

**Columns.** Each SAM field is a column of its own: FLAG, reference, POS, MAPQ,
CIGAR, the mate fields, SEQ, QUAL and the tags. ROOT splits each column into
pages of up to 1 MiB and compresses every page, with ZSTD by default, so values
of one kind are compressed together. Integer columns use RNTuple's split and
zigzag encodings, which put the bytes of equal significance next to each other.

**Mate fields.** PNEXT is stored as its distance from POS when the mate is on
the same reference, and TLEN as its difference from the value predicted from
the two positions, the strands and the CIGAR. For most pairs both are small
numbers that compress to almost nothing.

**Tags.** NM, AS, XS, MQ, RG, PG, MD, MC, SA and XA have typed columns of their
own; RG and PG values are stored once in METADATA and referred to by number.
Every other tag stays SAM text. Each record keeps the order its tags came in,
so they come back in that order.

**Blocks.** Records are grouped in blocks of 10,000 records or 5 million bases
or quality values, whichever comes first, as CRAM sizes its slices. The last
row of a block carries what is coded for the whole block:

- the read names, compressed with tok3, CRAM 3.1's name tokeniser;
- the quality scores, compressed with rANS order 1 as CRAM 3.1 compresses its
  quality series, or with fqzcomp (`-fqzcomp`);
- the consensus of the block's aligned reads, the commonest base at each
  position they cover. Each aligned read stores only its substitutions
  against the consensus and its soft-clipped and inserted bases, as CRAM's
  `embed_ref=2` slices store theirs;
- the block's reference: the consensus with the bases the reads' MD tags give
  where they differ from it, at a homozygous variant or a deletion, about 35
  bases a block of short reads. A read's MD and NM are computed again from it,
  as samtools calmd computes them, and stored only when that gives something
  else.

Unmapped reads and reads the consensus does not cover keep SEQ as text.

**METADATA.** The SAM header, the reference names, the records per reference
(for `ramstats idxstats`), where each block ends, the tag tables, whether the
file is sorted, and the format version.

**Region queries.** In a coordinate-sorted file, the first record of a region
is found by binary search over the positions, starting the longest reference
span in the file before the region; only the columns a query needs are read.

## Benchmarks

Two 5 Mb slices of chr20 (30,000,000 to 35,000,000 on GRCh38), one of short
reads and one of long reads:

| | Illumina | PacBio Revio |
|---|---|---|
| Sample | HG00154 from the 1000 Genomes 30x set (NYGC): NovaSeq, qualities binned to four values, Illumina read names | HG002 from PacBio's public Revio release, `HG002.m84011_220902_175841_s1.GRCh38.bam`: HiFi reads of about 15 kb |
| Records | 1,815,389 | 14,630 |
| Bases | 272 million | 228 million |

RAM is `bamtoramntuple -threads 4` with the default compression, 505 (zstd
level 10). CRAM is CRAM 3.1 with samtools' default profile, once against the
reference (`samtools view -O cram,version=3.1 -T chr20.fa`) and once with
`-O cram,version=3.1,embed_ref=2`, which stores a consensus of each slice's
reads instead, so that like RAM it reads without a reference. Run with
samtools 1.22 and ROOT 6.38 on an Intel i3-7020U laptop (2 cores, 4 threads).
[bench/README.md](bench/README.md) measures conversion, size and queries on
chr20 of NA12878 on your own machine.

### File size

| | Illumina | Revio |
|---|---:|---:|
| BAM | 152.4 MB | 101.2 MB |
| CRAM 3.1 with the reference | 67.9 MB | 70.1 MB |
| CRAM 3.1, `embed_ref=2` | 69.1 MB | 41.9 MB |
| RAM | 70.9 MB | 44.8 MB |
| RAM, `-compression 509` | 66.1 MB | 43.1 MB |

RAM gives both inputs back byte for byte (`ramdump -h` against
`samtools view -h`). Both CRAMs give the Illumina reads back exactly, but the
Revio reads with `M` in place of their `=` and `X` CIGAR operations and with
RG moved among the tags. Converted from the records as the CRAM returns them,
the Revio RAM file is 40.8 MB.

### By field

`ramstats size` for RAM and `samtools cram-size` for the `embed_ref=2` CRAM,
in MB on disk. CRAM keeps the CIGAR in its sequence features, so its SEQ row
holds both.

| | Illumina RAM | Illumina CRAM | Revio RAM | Revio CRAM |
|---|---:|---:|---:|---:|
| SEQ | 10.01 | 10.90 | 7.28 | 10.60 |
| CIGAR | 2.14 | | 5.66 | |
| QUAL | 11.05 | 11.07 | 24.54 | 24.51 |
| NAME | 7.38 | 7.38 | 0.05 | 0.06 |
| AUX (tags) | 34.06 | 33.94 | 7.20 | 6.52 |
| MATE (RNEXT, PNEXT, TLEN) | 3.87 | 2.78 | 0.03 | 0.03 |
| other (FLAG, POS, MAPQ, ...) | 2.31 | 2.71 | 0.03 | 0.06 |

Names and qualities are compressed with the codecs CRAM 3.1 uses, tok3 and
rANS order 1, in blocks the size of a CRAM slice, and come out the same size.
Neither stores MD and NM where it can compute them again: a CRAM from its
reference or consensus, RAM from each block's reference. On the Illumina reads
they took 5.0 MB when RAM stored them. Both store the bases as differences
from a consensus of the reads near them. On the Revio reads most of the CIGAR
is the `=` and `X` detail the CRAMs drop.

### Region queries

Records overlapping a region, counted (`samtools view -c`, `ramdump -c`) and
printed as SAM to `/dev/null`. CRAM is the `embed_ref=2` file. Best of five
runs, warm cache, in seconds; all three formats return the same records.

| Region | Records | Count BAM | Count CRAM | Count RAM | Print BAM | Print CRAM | Print RAM |
|---|---:|---:|---:|---:|---:|---:|---:|
| Illumina, 100 bp | 68 | 0.016 | 0.027 | 0.058 | 0.016 | 0.040 | 0.142 |
| Illumina, 85 kb | 23,408 | 0.030 | 0.041 | 0.059 | 0.044 | 0.098 | 0.201 |
| Illumina, 1 Mb | 726,726 | 0.840 | 0.698 | 0.114 | 1.700 | 3.050 | 1.983 |
| Revio, 100 bp | 20 | 0.012 | 0.017 | 0.041 | 0.012 | 0.042 | 0.107 |
| Revio, 85 kb | 180 | 0.022 | 0.022 | 0.041 | 0.024 | 0.070 | 0.110 |
| Revio, 1 Mb | 6,018 | 0.447 | 0.526 | 0.043 | 0.614 | 1.254 | 0.743 |

A RAM count reads only the position and CIGAR columns, so it stays fast as the
region grows. Printing a small region is slower than samtools: RAM decodes a
whole page of every column and a whole block of names, qualities and
consensus, up to 10,000 records or 5 million bases, about as much as a CRAM
slice.

### Compression settings

`-compression N` takes a ROOT compression code, algorithm times 100 plus level:
1 is ZLIB, 2 LZMA, 4 LZ4, 5 ZSTD, with levels 1 to 9, and 0 means no
compression. ROOT hands ZSTD twice the level, so the default 505 runs zstd at
level 10 and 509 at level 18. For a quick conversion use 501 or 503 with
`-threads`; for the smallest file, 509 at a much higher conversion cost. On the
two slices above, 509 makes the files 7% (Illumina) and 4% (Revio) smaller and
takes about four times as long to convert (56 s against 15 s for the Illumina
slice, 22 s against 5 s for Revio, `-threads 4`). Reading is no slower: a full
`ramdump` takes 3.7 s and 1.0 s with either.

`-fqzcomp` compresses the quality blocks with fqzcomp instead of rANS, CRAM
3.1's codec for its smaller profiles. On chr20 of NA12878 at 30x, whose
qualities are binned to four values, it saves about 1% of the qualities and
makes reading them three times as expensive: fqzcomp builds a 68 MB model for
every block it decodes. It may gain more on qualities that are not binned.
Counting records in a region reads no qualities and costs the same either way.

## Installation

RAM's tools build on Linux. The steps below are for Ubuntu; other
distributions need the same three things: ROOT, htslib and a C++17 compiler
with CMake 3.16 or newer.

### 1. System packages

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git pkg-config libhts-dev libtbb-dev libvdt-dev
```

`libhts-dev` is htslib, which reads the BAM input. `libtbb-dev` and `libvdt-dev`
are needed by ROOT. Read names and quality scores are compressed with
htscodecs, the codecs CRAM 3.1 uses; if `libhtscodecs-dev` is installed it is
used, otherwise the build downloads htscodecs 1.6.7 and compiles the files it
needs.

### 2. ROOT

RAM needs ROOT 6.38 or newer. The quickest way is the prebuilt archive from
the ROOT project; pick the one matching your Ubuntu version from
https://root.cern/install/all_releases/. For Ubuntu 24.04:

```bash
wget https://github.com/root-project/root/releases/download/v6-38-06/root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz
sudo tar -xzf root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz -C /opt/
```

ROOT has to be put on your path in every terminal you use it from:

```bash
source /opt/root/bin/thisroot.sh
```

Add that line to your `~/.bashrc` to make it permanent. `root --version` should
then print 6.38. If you use conda, `conda install -c conda-forge root` works as
well and needs no `thisroot.sh`.

### 3. Build

```bash
git clone https://github.com/Open-Genome/RAM.git
cd RAM
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The first configure downloads GoogleTest and Google Benchmark for the tests and
benchmarks, and htscodecs unless `libhtscodecs-dev` is installed. Without
network access, install `libhtscodecs-dev` and turn the tests and benchmarks
off with `-DRAMTOOLS_BUILD_TESTS=OFF -DRAMTOOLS_BUILD_BENCHMARKS=OFF`.

The tools are in `build/tools/`. Running one without arguments prints its
usage:

```bash
./build/tools/samtoramntuple
```

The tests take under a minute:

```bash
ctest --test-dir build
```

## Converting to RAM

From a SAM file:

```bash
./build/tools/samtoramntuple reads.sam reads.ram
```

From a BAM file:

```bash
./build/tools/bamtoramntuple reads.bam reads.ram
```

Without an output name the file is written next to the input, `reads.sam`
becoming `reads.ram`. Options for either tool:

| Option | Effect |
|--------|--------|
| `-compression N` | ROOT compression code, algorithm times 100 plus level; the default 505 is ZSTD, which ROOT runs at zstd level 10 (see [Compression settings](#compression-settings)) |
| `-threads N` | converts on N threads, see below |
| `-illumina` | stores quality scores in Illumina's 8 bins, which makes the file smaller |
| `-dropqual` | stores no quality scores |
| `-fqzcomp` | compresses quality scores with fqzcomp instead of rANS; not with `-illumina` or `-dropqual` |
| `-split` | `samtoramntuple` only: one file per reference, `reads_chr1.root`, `reads_chr2.root`, and so on |

`-illumina` and `-dropqual` lose information; without them a RAM file holds
everything the input did.

With `-threads N`, `samtoramntuple` reads the input in 64 MB blocks and every
thread parses, encodes and compresses blocks of its own. The blocks are written
in input order, so the records come out in the order they went in, whatever
the thread count. It needs about (N + 1) × 64 MB of memory on top of ROOT's
own: a block for every thread and one read ahead.
`bamtoramntuple` parses on one thread, codes each block of records (names,
qualities, the consensus) on a second and compresses on the others;
`samtoramntuple` with `-split` parses on one thread and compresses on the
others.

Sort the input by coordinate first, as `samtools sort` does. Region queries
find a region by binary search over the positions, which only works on a sorted
file. An unsorted file converts fine, but every query on it reads the whole
file; the converter warns when that happens.

## Querying a region

Regions are 1-based and inclusive, as in samtools:

| Region | Records |
|--------|---------|
| `chr1` | every record on `chr1` |
| `chr1:1000-2000` | records overlapping bases 1000 to 2000 |
| `chr1:1000` | records overlapping base 1000 only; samtools reads this as 1000 to the end of `chr1` |
| none, or `*` | every record in the file |

Reference names are whatever the input used, so a GRCh37 file is queried as
`13:32889611-32973805`, not `chr13`.

Count the records overlapping a region, as `samtools view -c` does:

```bash
./build/tools/ramdump -c reads.ram chr1:10150-10300
```

Print them as SAM lines:

```bash
./build/tools/ramdump reads.ram chr1:10150-10300
```

`ramdump` takes the `samtools view` options `-h` (include the header), `-H`
(header only), `-c` (count only), `-f` and `-F` (keep records with all, or none,
of these FLAG bits) and `-o FILE`. For example, primary alignments only:

```bash
./build/tools/ramdump -c -F 0x900 reads.ram chr1:10150-10300
```

`ramntupleview` counts as well and prints how long the query took:

```bash
./build/tools/ramntupleview reads.ram chr1:10150-10300
```

## Statistics

### Records per reference

`ramstats idxstats` prints what `samtools idxstats` prints: a line per
reference with its name, length, mapped and unmapped records, then `*` and the
records with no reference.

```bash
./build/tools/ramstats idxstats reads.ram
```

The converters store these counts in the file, so it reads no record and takes
the same time for any file. Unlike samtools, it does not need a sorted input.

### Records by FLAG

`ramstats flagstat` prints what `samtools flagstat` prints: the records counted
by their FLAG bits (primary, duplicates, mapped, properly paired, singletons and
so on), QC-passed and QC-failed apart. `-O tsv` and `-O json` give the other two
formats of samtools, which MultiQC reads.

```bash
./build/tools/ramstats flagstat reads.ram
```

It reads only the FLAG, reference, mate reference and MAPQ columns. On chr20 of
NA12878 at 30x (17.7 million records, 1000 Genomes high coverage), it takes
1.3 s of CPU on one thread. `samtools flagstat` takes 9.4 s on the BAM and
3.9 s on the CRAM 3.1, from which it also decodes only the fields it needs.

### Bytes per field

`ramstats size` prints how many bytes each field takes, the way
`samtools cram-size` does for a CRAM: a line per field with its bytes on disk,
its bytes before compression and its share of the file, then the total.

```bash
./build/tools/ramstats size reads.ram
```

It reads only the sizes RNTuple keeps for every page, so it takes the same
time for any file. `nameblock` and `qualblock` hold the read names and the
lossless qualities, already coded with tok3 and rANS (or fqzcomp) before
ROOT compresses them. SEQ is `consensus`, a block's consensus of its aligned
reads, with `record.seq_diff` and `record.seq_lit`, each read's substitutions
against it and its other bases. `consensus_patch` holds the reference bases MD
and NM are computed again from; a read's `record.tag_md` and `record.tag_nm`
are empty when they are.

## Checking a RAM file against samtools

A RAM file holds everything the input did, so it can be written back out and
compared:

```bash
./build/tools/ramdump -h reads.ram > roundtrip.sam
samtools view -c reads.bam chr1:10150-10300
./build/tools/ramdump -c reads.ram chr1:10150-10300
diff <(samtools idxstats reads.bam) <(./build/tools/ramstats idxstats reads.ram)
diff <(samtools flagstat reads.bam) <(./build/tools/ramstats flagstat reads.ram)
```

The two counts agree, the `diff`s print nothing, and for a SAM input
`roundtrip.sam` is the input byte for byte.

## Versioning

RAM follows [Semantic Versioning](https://semver.org). Until 1.0, a minor
release (0.2, 0.3, ...) may change the command-line options and the C++ API;
[CHANGELOG.md](CHANGELOG.md) lists what changed in each release.

Until 1.0, the file format may change between minor releases. A RAM file
records its format version; a reader opens only files of its own version and
tells you to convert again otherwise.

From 1.0 on, an option or public function that is going away is first marked
deprecated in the changelog and prints a warning for at least one minor release
before it is removed.

## Origin

RAM grew out of [ramtools](https://github.com/compiler-research/ramtools),
which began in 2017 at CERN openlab as a format stored in ROOT's TTree and was
rewritten on RNTuple as a Google Summer of Code project with CERN-HSF. This
repository continues that work on the format. To cite the ramtools
releases, use the "Cite this repository" button on the GitHub page, taken from
[CITATION.cff](CITATION.cff); every ramtools release is archived on Zenodo, and
[10.5281/zenodo.23362909](https://doi.org/10.5281/zenodo.23362909) stands for
all of them.

## Contributing

Bug reports, test cases and pull requests are welcome. [CONTRIBUTING.md](CONTRIBUTING.md)
describes the development build, the test suites, the formatting and clang-tidy
checks that run on every pull request, and what a reviewer looks for. The short
version: branch from `develop`, keep each PR to one topic, and put the samtools
command that confirms your result in the description.

## Reporting a problem

Open an issue at https://github.com/Open-Genome/RAM/issues. Include the
command you ran, its output, your ROOT version, and for a wrong result the
samtools command that gives the answer you expected. A few SAM lines that
reproduce the problem are the most useful thing you can attach.

## License

RAM is released under the Apache License 2.0; see [LICENSE](LICENSE).
