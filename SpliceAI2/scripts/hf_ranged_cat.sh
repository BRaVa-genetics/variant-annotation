#!/bin/bash
# usage: hf_ranged_cat.sh URL [PIECE_BYTES]   -> the file's bytes on stdout, fetched in byte ranges
# One long download of a 30 GB Hugging Face file can be cut short (on an HPC node, 2026-10-12: chr21 SNVs
# read 3.4M of 43.9M rows). Here each piece is a separate request on a fresh resolve (the signed CDN link
# expires), retried up to 6 times, written to local disk, size-checked, then emitted; the total must
# equal the file's size. Any failure exits non-zero, which pipefail turns into a failed task.
# The dataset is gated: accept its terms on Hugging Face, then put a read token in $HF_TOKEN or
# ~/.cache/huggingface/token.
set -euo pipefail
url=$1; step=${2:-1073741824}
tok=${HF_TOKEN:-$(tr -d '[:space:]' < ~/.cache/huggingface/token)}
size=$(curl -sSfIL -H "Authorization: Bearer $tok" "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{v=$2} END{print v+0}')
(( size > 0 )) || { echo "hf_ranged_cat: no size for $url" >&2; exit 1; }
tmp=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/hfcat.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
off=0; sent=0
while (( off < size )); do
  end=$(( off + step - 1 )); (( end > size - 1 )) && end=$(( size - 1 ))
  want=$(( end - off + 1 )); ok=""
  for try in 1 2 3 4 5 6; do
    if curl -sSfL -H "Authorization: Bearer $tok" -r "$off-$end" "$url" -o "$tmp/piece" \
       && [[ $(wc -c < "$tmp/piece") -eq $want ]]; then ok=1; break; fi
    echo "hf_ranged_cat: bytes $off-$end try $try failed, retrying" >&2; sleep $(( try * 10 ))
  done
  [[ -n $ok ]] || { echo "hf_ranged_cat: bytes $off-$end failed 6 times" >&2; exit 1; }
  cat "$tmp/piece"; sent=$(( sent + want )); off=$(( end + 1 ))
done
(( sent == size )) || { echo "hf_ranged_cat: sent $sent of $size bytes" >&2; exit 1; }
echo "hf_ranged_cat: $size bytes in $(( (size + step - 1) / step )) piece(s)" >&2
