# RAM

RAM (ROOT Alignment/Map) is a file format for aligned sequencing reads: the
same data as SAM or BAM, stored column by column in ROOT's RNTuple. A RAM file
gives back every record byte for byte, needs no reference to be read, answers
region queries, and is less than half the size of the BAM:

| 5 Mb of chr20 | BAM | CRAM 3.1, `embed_ref=2` | RAM |
|---|---:|---:|---:|
| Illumina short reads (HG00154, 30x) | 152.4 MB | 69.1 MB | 70.9 MB |
| PacBio Revio HiFi reads (HG002) | 101.2 MB | 41.9 MB | 44.8 MB |

| Tool | What it does |
|------|--------------|
| `samtoramntuple`, `bamtoramntuple` | convert SAM or BAM to RAM |
| `ramdump` | writes a RAM file, or a region of it, as SAM; takes the `samtools view` options |
| `ramstats` | `idxstats`, `flagstat` (as samtools prints them) and `size`, the bytes each field takes |
| `ramntupleview` | counts the records in a region and times the query |

## Quick start

RAM's tools need ROOT 6.38 or newer, htslib and a C++17 compiler with CMake
3.16 or newer. On Ubuntu 24.04:

```bash
sudo apt-get install -y build-essential cmake git pkg-config libhts-dev libtbb-dev libvdt-dev
wget https://github.com/root-project/root/releases/download/v6-38-06/root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz
sudo tar -xzf root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz -C /opt/
source /opt/root/bin/thisroot.sh        # in every terminal, or in ~/.bashrc

git clone https://github.com/Open-Genome/RAM.git
cd RAM
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
ctest --test-dir build                  # under a minute
```

Other Ubuntu versions have their own ROOT archive at
https://root.cern/install/all_releases/; with conda,
`conda install -c conda-forge root` replaces the ROOT steps. The first
configure downloads GoogleTest, Google Benchmark and htscodecs (unless
`libhtscodecs-dev` is installed); offline, install `libhtscodecs-dev` and add
`-DRAMTOOLS_BUILD_TESTS=OFF -DRAMTOOLS_BUILD_BENCHMARKS=OFF`.

Convert a coordinate-sorted BAM, query a region and check the result against
samtools:

```bash
./build/tools/bamtoramntuple reads.bam reads.ram -threads 4
./build/tools/ramdump reads.ram chr1:10150-10300
./build/tools/ramdump -c reads.ram chr1:10150-10300        # = samtools view -c
diff <(./build/tools/ramdump -h reads.ram) <(samtools view -h --no-PG reads.bam)
```

## How RAM stores reads

- **Columns.** Every SAM field is a column; ROOT compresses it in pages of up
  to 1 MiB, with ZSTD by default. A mate's PNEXT on the same reference is
  stored as its distance from POS, and TLEN as its difference from the value
  predicted from the pair's positions, so for most pairs both are near zero.
- **Tags.** NM, AS, XS, MQ, RG, PG, MD, MC, SA and XA have typed columns;
  other tags stay SAM text. Tags come back in their input order.
- **Blocks.** Every 10,000 records, or 5 million bases, form a block, as a
  CRAM slice does. A block stores its read names with tok3 and its qualities
  with rANS order 1, the codecs of CRAM 3.1, and the consensus of its reads:
  each aligned read keeps only its substitutions against it and its clipped
  and inserted bases. MD and NM are computed again against the block's
  reference (the consensus, corrected at the few bases where the MD tags show
  it differs) and stored only when that gives something else.
- **Metadata.** The SAM header, reference names, records per reference, block
  boundaries and the format version.

Region queries find their first record by binary search over the positions,
so convert coordinate-sorted input: an unsorted file converts, but every query
reads it in full.

## Benchmarks

Two 5 Mb slices of chr20 (30–35 Mb on GRCh38):

| | Illumina | PacBio Revio |
|---|---|---|
| Sample | HG00154, 1000 Genomes 30x (NYGC): NovaSeq, four quality bins | HG002, PacBio's public Revio release: HiFi reads of about 15 kb |
| Records | 1,815,389 | 14,630 |
| Bases | 272 million | 228 million |

RAM is `bamtoramntuple -threads 4` at the default compression. CRAM is CRAM
3.1 from samtools 1.22, with the reference and with `embed_ref=2`, which
embeds a consensus so that, like RAM, it reads without the reference. Intel
i3-7020U laptop (2 cores, 4 threads), ROOT 6.38.
[bench/README.md](bench/README.md) runs a benchmark on your own machine.

### Size

MB on disk, per field from `ramstats size` and `samtools cram-size`. CRAM
counts the CIGAR under SEQ.

| | Illumina RAM | Illumina CRAM | Revio RAM | Revio CRAM |
|---|---:|---:|---:|---:|
| SEQ | 10.01 | 10.90 | 7.28 | 10.60 |
| CIGAR | 2.14 | | 5.66 | |
| QUAL | 11.05 | 11.07 | 24.54 | 24.51 |
| NAME | 7.38 | 7.38 | 0.05 | 0.06 |
| AUX (tags) | 34.06 | 33.94 | 7.20 | 6.52 |
| MATE (RNEXT, PNEXT, TLEN) | 3.87 | 2.78 | 0.03 | 0.03 |
| other (FLAG, POS, MAPQ, ...) | 2.31 | 2.71 | 0.03 | 0.06 |
| **File size** | **70.9** | **69.1** | **44.8** | **41.9** |

CRAM with the reference is 67.9 MB (Illumina) and 70.1 MB (Revio); BAM is
152.4 and 101.2 MB. RAM gives both inputs back byte for byte; both CRAMs give
the Revio reads back with `M` in place of their `=` and `X` CIGAR operations,
most of RAM's 5.66 MB of Revio CIGAR. Converted from the records the CRAM
returns, the Revio RAM file is 40.8 MB.

### Region queries

Seconds, best of five, warm cache; CRAM is the `embed_ref=2` file. All three
return the same records.

| Region | Records | Count BAM | Count CRAM | Count RAM | Print BAM | Print CRAM | Print RAM |
|---|---:|---:|---:|---:|---:|---:|---:|
| Illumina, 100 bp | 68 | 0.016 | 0.027 | 0.058 | 0.016 | 0.040 | 0.142 |
| Illumina, 85 kb | 23,408 | 0.030 | 0.041 | 0.059 | 0.044 | 0.098 | 0.201 |
| Illumina, 1 Mb | 726,726 | 0.840 | 0.698 | 0.114 | 1.700 | 3.050 | 1.983 |
| Revio, 100 bp | 20 | 0.012 | 0.017 | 0.041 | 0.012 | 0.042 | 0.107 |
| Revio, 85 kb | 180 | 0.022 | 0.022 | 0.041 | 0.024 | 0.070 | 0.110 |
| Revio, 1 Mb | 6,018 | 0.447 | 0.526 | 0.043 | 0.614 | 1.254 | 0.743 |

A count reads only the position and CIGAR columns, so it stays fast as the
region grows. Printing a small region is slower than samtools because RAM
decodes a whole block of names, qualities and consensus, about a CRAM slice.

## Converting

```bash
./build/tools/samtoramntuple reads.sam [reads.ram] [options]
./build/tools/bamtoramntuple reads.bam [reads.ram] [options]
```

Without an output name, `reads.sam` or `reads.bam` becomes `reads.ram`.

| Option | Effect |
|--------|--------|
| `-threads N` | converts on N threads |
| `-compression N` | ROOT's code, algorithm × 100 + level (1 ZLIB, 2 LZMA, 4 LZ4, 5 ZSTD; 0 none). The default 505 runs zstd at level 10, since ROOT doubles the level |
| `-fqzcomp` | compresses qualities with fqzcomp instead of rANS |
| `-illumina` | bins qualities to Illumina's 8 values (lossy) |
| `-dropqual` | stores no qualities (lossy) |
| `-split` | `samtoramntuple` only: one file per reference, `reads_chr1.root` and so on |

`-compression 509` makes the benchmark files 7% (Illumina) and 4% (Revio)
smaller, takes about four times as long to convert, and reads as fast.
`-fqzcomp` saved about 1% of the qualities of binned NA12878 reads and makes
reading them three times as expensive; it may gain more on unbinned
qualities. With `-threads N`, `samtoramntuple` parses in 64 MB blocks on every
thread and needs about (N + 1) × 64 MB more memory; the records keep their
input order.

## Querying

Regions are 1-based and inclusive, as in samtools; reference names are the
input's (`13:32889611-32973805` in a GRCh37 file).

| Region | Records |
|--------|---------|
| `chr1` | every record on `chr1` |
| `chr1:1000-2000` | records overlapping bases 1000 to 2000 |
| `chr1:1000` | records overlapping base 1000 (samtools reads this as 1000 to the end) |
| none, or `*` | every record |

`ramdump` takes the `samtools view` options `-h`, `-H`, `-c`, `-f`, `-F` and
`-o FILE`:

```bash
./build/tools/ramdump -c -F 0x900 reads.ram chr1:10150-10300   # primary alignments only
```

`ramstats` prints statistics without reading records, or only the columns it
needs:

```bash
./build/tools/ramstats idxstats reads.ram     # = samtools idxstats, from the counts stored in the file
./build/tools/ramstats flagstat reads.ram     # = samtools flagstat; -O tsv and -O json as well
./build/tools/ramstats size reads.ram         # bytes per field, like samtools cram-size
```

In `ramstats size`, `nameblock` and `qualblock` are the names and qualities;
SEQ is `consensus`, `record.seq_diff` and `record.seq_lit`; `consensus_patch`
holds the reference bases MD and NM are computed against.

## Versioning

RAM follows [Semantic Versioning](https://semver.org); [CHANGELOG.md](CHANGELOG.md)
lists every change. Until 1.0 the options, the C++ API and the file format may
change between minor releases: a file records its format version, and a
reader of another version asks you to convert it again.

## Origin, contributing and license

RAM grew out of [ramtools](https://github.com/compiler-research/ramtools),
which began in 2017 at CERN openlab on ROOT's TTree and was rewritten on
RNTuple as a Google Summer of Code project with CERN-HSF. Its releases are
archived on Zenodo under [10.5281/zenodo.23362909](https://doi.org/10.5281/zenodo.23362909);
[CITATION.cff](CITATION.cff) gives the reference.

Issues and pull requests are welcome at https://github.com/Open-Genome/RAM;
[CONTRIBUTING.md](CONTRIBUTING.md) describes the checks every pull request
runs. For a wrong result, include the command, its output and the samtools
command that gives the answer you expected.

RAM is released under the Apache License 2.0; see [LICENSE](LICENSE).
