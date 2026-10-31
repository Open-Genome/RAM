#!/bin/sh
# Writes RAM files with ROOT and checks that ramread prints exactly what ramdump
# prints, as SAM, as BAM read back by samtools, and for a region, so a change to
# the format or to ROOT's encodings fails here first.
# Usage: check.sh <directory with samtoramntuple and ramdump> <ramread>
set -eu
tools=$1
ramread=$2
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT

# Paired reads on two references, with unmapped reads, mates on the other
# reference, missing qualities and array tags.
awk 'BEGIN {
   srand(7)
   printf "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrA\tLN:5000000\n@SQ\tSN:chrB\tLN:5000000\n@RG\tID:g\tSM:s\n"
   split("A C G T", base, " ")
   for (i = 0; i < 40000; ++i) {
      if (i % 20000 == 0) pos = 1
      pos += int(rand() * 200)
      seq = qual = ""
      for (j = 0; j < 150; ++j) {
         seq = seq base[int(rand() * 4) + 1]
         qual = qual sprintf("%c", 35 + int(rand() * 7) * 5)
      }
      ref = i < 20000 ? "chrA" : "chrB"
      if (i % 97 == 0)
         printf "r%d\t69\t%s\t%d\t0\t*\t=\t%d\t0\t%s\t%s\n", i, ref, pos, pos, seq, qual
      else if (i % 89 == 0)
         printf "r%d\t65\t%s\t%d\t60\t150M\t%s\t%d\t0\t%s\t*\tNM:i:0\n", i, ref, pos, ref == "chrA" ? "chrB" : "chrA", pos, seq
      else
         printf "r%d\t99\t%s\t%d\t60\t100M2I48M\t=\t%d\t%d\t%s\t%s\tNM:i:2\tMD:Z:148\tAS:i:%d\tRG:Z:g\tZB:B:i,-3,%d\n",
            i, ref, pos, pos + 250, 400, seq, qual, int(rand() * 150), i
   }
}' > "$dir/in.sam"

for opts in "" -fqzcomp -illumina -dropqual "-threads 4" -split; do
   rm -f "$dir"/*.ram "$dir"/*.root
   # shellcheck disable=SC2086
   "$tools/samtoramntuple" "$dir/in.sam" "$dir/out.ram" $opts > /dev/null
   for f in "$dir"/*.ram "$dir"/*.root; do
      [ -e "$f" ] || continue
      "$tools/ramdump" -h "$f" > "$dir/want.sam"
      "$ramread" -h "$f" > "$dir/got.sam"
      cmp "$dir/want.sam" "$dir/got.sam" || { echo "ramread differs from ramdump: $opts ${f##*/}"; exit 1; }
      "$ramread" -b "$f" | samtools view -h --no-PG - > "$dir/got.sam"
      cmp "$dir/want.sam" "$dir/got.sam" || { echo "ramread -b differs from ramdump: $opts ${f##*/}"; exit 1; }
      "$tools/ramdump" -h "$f" chrB:1000000-1500000 > "$dir/want.sam"
      "$ramread" -h "$f" chrB:1000000-1500000 > "$dir/got.sam"
      cmp "$dir/want.sam" "$dir/got.sam" || { echo "ramread differs from ramdump on a region: $opts ${f##*/}"; exit 1; }
   done
   echo "ok: ${opts:-default}"
done
