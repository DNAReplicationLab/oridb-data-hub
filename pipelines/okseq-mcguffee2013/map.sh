#!/bin/bash
# Ok-seq mapping: adapter/quality trim (fastp) -> bowtie2 (sacCer3 + 2micron decoy)
#                 -> proper pairs only, MAPQ>=20, read1 records kept (TLEN gives the fragment)
set -eu -o pipefail
R1=$1; R2=$2; OUT=$3; BT=${4:-9}
IDX="$WS/okseq/ref/sacCer3_plus2mu"
fastp -i "$R1" -I "$R2" --stdout --detect_adapter_for_pe \
      -l 20 -q 15 -u 40 -w 3 -j "${OUT}.fastp.json" -h /dev/null 2> "${OUT}.fastp.log" \
| bowtie2 -p $BT -x "$IDX" --interleaved - --no-unal -X 2000 2> "${OUT}.bowtie2.log" \
| samtools view -@2 -b -f 67 -F 3852 -q 20 - \
| samtools sort -@2 -m 1500M -o "${OUT}.r1.bam" -
samtools index "${OUT}.r1.bam"
samtools idxstats "${OUT}.r1.bam" > "${OUT}.idxstats.txt"
echo "DONE $OUT"
