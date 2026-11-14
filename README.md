# RAM

RAM (ROOT Alignment/Map) is a file format for aligned sequencing reads. It
holds the same data as SAM or BAM, stored column by column in ROOT's RNTuple:
every record comes back byte for byte, no reference is needed to read it, and
regions can be queried as with samtools.

This repository has the tools to convert SAM and BAM files to RAM, query them
and write them back out as SAM.

## Install

You need ROOT 6.38 or newer, htslib, CMake 3.16 or newer and a C++17
compiler. The steps below are for Ubuntu 24.04.

**1. Install the system packages**

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git pkg-config libhts-dev libtbb-dev libvdt-dev
```

**2. Install ROOT**

```bash
wget https://github.com/root-project/root/releases/download/v6-38-06/root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz
sudo tar -xzf root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz -C /opt/
source /opt/root/bin/thisroot.sh
```

Run `source /opt/root/bin/thisroot.sh` in every new terminal, or add it to
`~/.bashrc`. For other Ubuntu versions pick the matching archive from
https://root.cern/install/all_releases/; with conda,
`conda install -c conda-forge root` replaces this step.

**3. Build RAM**

```bash
git clone https://github.com/Open-Genome/RAM.git
cd RAM
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

**4. Run the tests** (under a minute)

```bash
ctest --test-dir build
```

The tools are now in `build/tools/`. Run any of them without arguments to see
its usage.

The first build downloads GoogleTest, Google Benchmark and htscodecs. Without
network access, install `libhtscodecs-dev` and add
`-DRAMTOOLS_BUILD_TESTS=OFF -DRAMTOOLS_BUILD_BENCHMARKS=OFF` to the first
`cmake` command.

## Convert to RAM

Sort the input by coordinate first, so that region queries are fast:

```bash
samtools sort -o reads.sorted.bam reads.bam
```

Then convert it:

```bash
./build/tools/bamtoramntuple reads.sorted.bam reads.ram -threads 4    # from BAM
./build/tools/samtoramntuple reads.sam reads.ram -threads 4           # from SAM
```

The input comes first, then the output name (optional: `reads.bam` becomes
`reads.ram`), then the options:

| Option | What it does |
|--------|--------------|
| `-threads N` | converts on N threads |
| `-compression N` | ROOT compression code; the default `505` is ZSTD, `509` gives a smaller file but converts about four times slower |
| `-fqzcomp` | compresses quality scores with fqzcomp instead of rANS |
| `-illumina` | bins quality scores to Illumina's 8 values (loses information) |
| `-dropqual` | drops quality scores (loses information) |
| `-split` | `samtoramntuple` only: one file per reference (`reads_chr1.root`, ...) |

## Query a region

Regions are written as in samtools and are 1-based and inclusive. Use the
reference names of your input (`chr1` or `1`).

```bash
./build/tools/ramdump reads.ram chr1:10000-20000           # print the records as SAM
./build/tools/ramdump -c reads.ram chr1:10000-20000        # count them
./build/tools/ramdump -h reads.ram > reads.sam             # whole file, with the header
./build/tools/ramdump -H reads.ram                         # header only
./build/tools/ramdump -F 0x900 reads.ram chr1              # primary alignments on chr1
./build/tools/ramdump -o out.sam reads.ram chr1:10000-20000   # write to a file
./build/tools/ramntupleview reads.ram chr1:10000-20000     # count and time the query
```

`ramdump` takes the `samtools view` options `-h`, `-H`, `-c`, `-f`, `-F` and
`-o`. One difference: `chr1:1000` means base 1000 only, where samtools reads it
as 1000 to the end of `chr1`.

## Statistics

```bash
./build/tools/ramstats idxstats reads.ram        # records per reference, as samtools idxstats
./build/tools/ramstats flagstat reads.ram        # as samtools flagstat
./build/tools/ramstats flagstat -O json reads.ram   # also -O tsv
./build/tools/ramstats size reads.ram            # bytes each field takes, as samtools cram-size
```

## Check a RAM file against samtools

Every result can be compared with samtools. None of these print anything when
the RAM file matches:

```bash
diff <(./build/tools/ramdump -h reads.ram) <(samtools view -h --no-PG reads.sorted.bam)
diff <(./build/tools/ramstats idxstats reads.ram) <(samtools idxstats reads.sorted.bam)
diff <(./build/tools/ramstats flagstat reads.ram) <(samtools flagstat reads.sorted.bam)
diff <(./build/tools/ramdump -c reads.ram chr1:10000-20000) <(samtools view -c reads.sorted.bam chr1:10000-20000)
```

`samtools idxstats` and region queries need the BAM index:
`samtools index reads.sorted.bam`.

## Versioning

RAM follows [Semantic Versioning](https://semver.org); [CHANGELOG.md](CHANGELOG.md)
lists every change. Until 1.0 the file format may change between releases: a
RAM file records its format version, and a newer `ramdump` asks you to convert
an older file again.

## Origin, contributing and license

RAM grew out of [ramtools](https://github.com/compiler-research/ramtools),
which began in 2017 at CERN openlab and was rewritten on RNTuple as a Google
Summer of Code project with CERN-HSF. Its releases are archived on Zenodo
under [10.5281/zenodo.23362909](https://doi.org/10.5281/zenodo.23362909).

Issues and pull requests are welcome at https://github.com/Open-Genome/RAM;
see [CONTRIBUTING.md](CONTRIBUTING.md).

RAM is released under the Apache License 2.0; see [LICENSE](LICENSE).
