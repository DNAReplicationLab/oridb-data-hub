#!/bin/bash
# Wait for each run's chunk set to be complete, assemble, verify md5, map.
set -u
cd "$WS/okseq/fastq"
LOG="$WS/okseq/pipeline.log"; : > "$LOG"
complete() {  # $1=name $2=expected total bytes
  local name=$1 exp=$2 tot
  tot=$(cat parts/$name/* 2>/dev/null | wc -c)
  [ "$tot" = "$exp" ]
}
sizeof() { cat parts/$1/* 2>/dev/null | wc -c; }
run_one() {
  local run=$1 l1=$2 l2=$3 md1=$4 md2=$5
  while true; do
    s1=$(du -sk parts/${run}_1 | awk '{print $1}'); s2=$(du -sk parts/${run}_2 | awk '{print $1}')
    echo "$(date +%H:%M:%S) $run waiting: _1 $((s1/1024))MB/$((l1/1048576))MB  _2 $((s2/1024))MB/$((l2/1048576))MB" >> "$LOG"
    if [ $((s1*1024)) -ge $((l1-2000000)) ] && [ $((s2*1024)) -ge $((l2-2000000)) ]; then
      sleep 20
      cat parts/${run}_1/* > ${run}_1.asm.fastq.gz; cat parts/${run}_2/* > ${run}_2.asm.fastq.gz
      a1=$(wc -c < ${run}_1.asm.fastq.gz); a2=$(wc -c < ${run}_2.asm.fastq.gz)
      if [ "$a1" = "$l1" ] && [ "$a2" = "$l2" ]; then
        o1=$(md5 -q ${run}_1.asm.fastq.gz); o2=$(md5 -q ${run}_2.asm.fastq.gz)
        echo "$(date +%H:%M:%S) $run assembled. md5_1 $( [ "$o1" = "$md1" ] && echo OK || echo MISMATCH:$o1 ) md5_2 $( [ "$o2" = "$md2" ] && echo OK || echo MISMATCH:$o2 )" >> "$LOG"
        break
      fi
      echo "$(date +%H:%M:%S) $run size mismatch ($a1/$l1 $a2/$l2), still downloading" >> "$LOG"
    fi
    sleep 60
  done
  echo "$(date +%H:%M:%S) $run mapping start" >> "$LOG"
  bash "$WS/okseq/map.sh" ${run}_1.asm.fastq.gz ${run}_2.asm.fastq.gz "$WS/okseq/bam/$run" 9 >> "$LOG" 2>&1
  echo "$(date +%H:%M:%S) $run mapping done" >> "$LOG"
}
mkdir -p "$WS/okseq/bam"
run_one SRR566701 8655237045 9051017187 ddfef82d1978526f99535ce5fc314940 04b93503b9060314c663fb4934cc25b6
run_one SRR566702 4260407357 4240339756 d05f3bd5c43903f0497265bcda1389da 583f92fbffc232ea6f71569f7f255ef0
echo "ALL DONE" >> "$LOG"
tail -5 "$LOG"
