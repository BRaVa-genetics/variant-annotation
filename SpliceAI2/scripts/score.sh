#!/bin/bash
# usage: [NPROC=4] score.sh WORK_DIR CHUNK.tsv CONTROLS.tsv      (on a GPU node; WORK_DIR from setup.sh)
# Scores CHUNK.tsv (gpu_set.py: chrom pos ref alt strand, header) with Illumina's SpliceAI2 scorer and writes
# CHUNK.scores.tsv.gz (chrom pos ref alt strand spliceai2_summary_score), which brava_create_annot.py --spliceai2
# reads next to the lookup output. The 300 control variants (control.py pick) are scored in the same run, and
# nothing is written unless they agree with Illumina's precomputed scores (control.py compare): a bad node or GPU
# cannot leave scores behind. A finished chunk is skipped.
# Speed: the scorer runs one variant at a time with its main Python thread at 100% of a core (A100: ~4.4 variants/s
# per process), so NPROC copies run side by side on the one GPU and their outputs are joined in input order;
# NPROC=4 on an A100: 1,278 variants in 7 min 20 s, model loading included (2026-10-17), with summary scores
# identical to those of the faster in-house scorer used for BRaVa's October 2026 release, on all 978 non-control variants.
set -euo pipefail
W=$(cd "$1" && pwd); in=$2; ctl=$3; out=${in%.tsv}.scores.tsv.gz; NPROC=${NPROC:-4}
here=$(cd "$(dirname "$0")" && pwd)
[[ -s $out ]] && { echo "$in: done already"; exit 0; }
command -v nvidia-smi >/dev/null || { echo "no GPU on $(hostname): run inside a GPU allocation"; exit 1; }
ENV=$W/sai2env; export HF_HOME=$W/.cache/hf TORCHINDUCTOR_CACHE_DIR=$W/.cache/inductor PYTHONFAULTHANDLER=1
work=${in%.tsv}.work.$$; mkdir -p "$work"
n=$(( $(wc -l < "$in") - 1 )); m=$(( $(wc -l < "$ctl") - 1 ))
{ cat "$in"; tail -n +2 "$ctl" | cut -f1-5; } > "$work/all.tsv"
per=$(( (n + m + NPROC - 1) / NPROC ))
tail -n +2 "$work/all.tsv" | split -l "$per" -d -a 2 - "$work/p"
pids=()
for f in "$work"/p??; do
  { head -1 "$in"; cat "$f"; } > "$f.tsv"; rm "$f"
  "$ENV/bin/spliceai2" --model_folders "$W"/models/SpliceAI2/spliceai2_multispecies_13M_01 "$W"/models/SpliceAI2/spliceai2_multispecies_13M_02 \
      --var_tsv_file "$f.tsv" --fasta_file "$W/ref/hg38.fa" > "$f.log" 2>&1 &
  pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
(( fail == 0 )) || { echo "a spliceai2 process failed:"; tail -5 "$work"/p??.tsv.log; exit 1; }
# the scorer writes 260 columns per variant; keep the key and the summary score (max of the four site deltas)
first=1
for f in "$work"/p??.tsv; do
  [[ -s $f.spliceai2 ]] || { echo "no output $f.spliceai2"; exit 1; }
  awk -F'\t' -v first=$first 'NR==1{for(i=1;i<=NF;i++) ix[$i]=i; if(first) print "chrom\tpos\tref\talt\tstrand\tspliceai2_summary_score"; next}
    {print $ix["chrom"]"\t"$ix["pos"]"\t"$ix["ref"]"\t"$ix["alt"]"\t"$ix["strand"]"\t"$ix["spliceai2_summary_score"]}' "$f.spliceai2"
  first=0
done > "$work/scored.tsv"
k=$(( $(wc -l < "$work/scored.tsv") - 1 ))
[[ $k -eq $(( n + m )) ]] || { echo "$in: $k rows scored, expected $(( n + m )) -- not accepted"; exit 1; }
{ head -1 "$work/scored.tsv"; tail -n "$m" "$work/scored.tsv"; } | gzip > "$work/ctl.scores.tsv.gz"
python3 "$here/control.py" compare "$ctl" "$work/ctl.scores.tsv.gz" || { echo "$in: controls FAILED -- not accepted"; exit 1; }
head -n $(( n + 1 )) "$work/scored.tsv" | gzip > "$out.tmp" && mv "$out.tmp" "$out" && rm -rf "$work"
echo "$in: $n variants scored, controls PASS -> $out"
