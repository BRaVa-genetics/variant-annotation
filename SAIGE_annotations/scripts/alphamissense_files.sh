#!/bin/bash
# usage: alphamissense_files.sh UNIVERSAL_SAIGE_DIR GROUP_FILE.txt.gz CHROM [SPLICE_LABEL=spliceai2_gt040]
# The two AlphaMissense files of one chromosome, from brava_create_annot.py's group file and long table
# (GROUP_FILE.long.tsv.gz):
#   GROUP.flexrv_AMPCT.txt.gz         the group file plus a 'score:AMPCT' line per gene, for FlexRV (universal-saige
#                                     step 2 --flexRVscore AMPCT), on the FlexRV paper's scale (Schwartzentruber et al.
#                                     2025): pLoF 1.0; a missense variant the percentile of its AlphaMissense score among
#                                     all possible missense variants of the gene's transcript, (#below + (#equal+1)/2)/N;
#                                     protein-altering non-missense entries (in-frame indels, LOFTEE-LC LoF, start/stop
#                                     lost) the gene's mean, after ranking; entries with no protein change dropped
#   GROUP.am_moderate.saige_group.txt.gz   missense variants with raw AlphaMissense >= 0.906 (ACMG/AMP PP3 Moderate),
#                                     label am_moderate, for a separate SKAT-O run
# Both come from ONE run of universal-saige's 04_flexrv_groupfile.sh (the canonical FlexRV builder), each entry scored on
# its own transcript (--annoTable: MANE Select, else canonical, the version GENCODE v39 holds): the am_moderate mask reads
# the raw AlphaMissense score of the row each entry's percentile came from (since 2026-10-23; the 2026-10 release's masks
# came from a separate raw run). chrX also reads AlphaMissense's chrY PAR rows, and paralog copies AlphaMissense does not
# score take its protein-level scores when universal-saige's alphamissense_paralogs.py output sits next to the
# per-chromosome files. KEEP_RAW=1 also writes flexrv_AM.txt.gz, the raw-score FlexRV weights of the earlier release.
# The full mapping from raw AlphaMissense to the percentile is in universal-saige's README, step 4. Its AlphaMissense,
# MANE and GENCODE v39 files: bash download_resources.sh --alphamissense in universal-saige. Every missense-type label
# can hold a missense variant -- including the SpliceAI category (a missense variant with low REVEL/CADD and a splice
# call) -- so all three are passed to --missenseAnno, and the FlexRV mask in step 2 is then
# pLoF:damaging_missense:other_missense:SPLICE_LABEL.
set -euo pipefail
US=$(cd "$1" && pwd); g=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); c=${3#chr}; lab=${4:-spliceai2_gt040}
here=$(cd "$(dirname "$0")" && pwd)
stem=${g%.txt.gz}; stem=${stem%.gz}
R=$US/resources
am=$R/AlphaMissense_hg38.tsv.gz; iso=$R/AlphaMissense_isoforms_hg38.tsv.gz
[[ -s $R/am_by_chrom/AlphaMissense_hg38.chr$c.tsv.gz ]] && am=$R/am_by_chrom/AlphaMissense_hg38.chr$c.tsv.gz
[[ -s $R/am_by_chrom/AlphaMissense_isoforms_hg38.chr$c.tsv.gz ]] && iso=$R/am_by_chrom/AlphaMissense_isoforms_hg38.chr$c.tsv.gz
[[ -s $am && -s $R/gencode.v39.annotation.gtf.gz && -s $R/MANE.GRCh38.summary.txt.gz ]] \
  || { echo "AlphaMissense / GENCODE v39 / MANE files missing in $R: run 'bash download_resources.sh --alphamissense' in universal-saige"; exit 1; }
[[ -s $g.long.tsv.gz ]] || { echo "no long table $g.long.tsv.gz next to the group file"; exit 1; }
common=(--group "$g" --chr "$c" --missenseAnno damaging_missense,other_missense,$lab
        --am "$am" --isoforms "$iso" --annoTable "$g.long.tsv.gz" --mane "$R/MANE.GRCh38.summary.txt.gz"
        --gtf "$R/gencode.v39.annotation.gtf.gz")
# the FlexRV paper's scale: the FlexRV weights, and the provenance the am_moderate mask is read from
( cd "$US" && bash 04_flexrv_groupfile.sh "${common[@]}" --name AMPCT --amScale percentile --nonMissense gene-mean \
    --provenance "$stem.flexrv_AMPCT.provenance.tsv" --out "$stem.flexrv_AMPCT.txt" )
gzip -f "$stem.flexrv_AMPCT.txt" "$stem.flexrv_AMPCT.provenance.tsv"
python3 "$here/am_mask.py" "$g" "$stem.flexrv_AMPCT.provenance.tsv.gz" "$stem.am_moderate.saige_group.txt.gz"
if [[ ${KEEP_RAW:-0} == 1 ]]; then
  ( cd "$US" && bash 04_flexrv_groupfile.sh "${common[@]}" --name AM --nonMissense drop --out "$stem.flexrv_AM.txt" )
  gzip -f "$stem.flexrv_AM.txt"
fi
echo "wrote $stem.flexrv_AMPCT.txt.gz and $stem.am_moderate.saige_group.txt.gz"
