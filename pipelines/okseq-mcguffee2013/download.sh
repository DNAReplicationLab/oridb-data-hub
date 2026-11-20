#!/bin/bash
# Resumable chunked fetch from ENA. Each chunk is a byte range appended to until complete;
# stalled connections are aborted by --speed-limit/--speed-time and retried.
set -u
cd "$WS/okseq/fastq"
BASE=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566
CHUNK=200000000
rm -f jobs.txt fails.txt
while read -r run name len; do
  mkdir -p parts/$name
  n=$(( (len + CHUNK - 1) / CHUNK ))
  for ((i=0;i<n;i++)); do
    s=$((i*CHUNK)); e=$((s+CHUNK-1)); [ $e -ge $len ] && e=$((len-1))
    printf "%s\t%s\t%s\t%s\n" "$BASE/$run/$name.fastq.gz" "parts/$name/$(printf %04d $i)" "$s" "$e" >> jobs.txt
  done
done <<'TBL'
SRR566701 SRR566701_1 8655237045
SRR566701 SRR566701_2 9051017187
SRR566702 SRR566702_1 4260407357
SRR566702 SRR566702_2 4240339756
TBL
echo "chunks: $(wc -l < jobs.txt)"
fetch() {
  f=$1; o=$2; s=$3; e=$4; exp=$((e-s+1))
  for t in $(seq 1 40); do
    got=0; [ -f "$o" ] && got=$(wc -c < "$o")
    if [ "$got" -ge "$exp" ]; then [ "$got" -gt "$exp" ] && : > "$o" && continue; exit 0; fi
    curl -s --fail --speed-limit 30000 --speed-time 45 --max-time 1200 \
         -r $((s+got))-$e "$f" >> "$o"
    sleep 2
  done
  echo "FAIL $o" >> fails.txt
}
export -f fetch
cat jobs.txt | xargs -P 12 -n 4 bash -c 'fetch "$1" "$2" "$3" "$4"' _
for name in SRR566701_1 SRR566701_2 SRR566702_1 SRR566702_2; do
  cat parts/$name/* > $name.fastq.gz && echo "$name assembled $(wc -c < $name.fastq.gz)"
done
[ -f fails.txt ] && cat fails.txt || echo "no chunk failures"
md5 -q SRR566701_1.fastq.gz SRR566701_2.fastq.gz SRR566702_1.fastq.gz SRR566702_2.fastq.gz > md5_observed.txt
paste <(cut -d' ' -f1 ena_md5.txt) md5_observed.txt | awk '{print ($1==$2?"MD5_OK":"MD5_MISMATCH"), $1, $2}'
