#!/bin/bash
# usage: alphamissense_files.sh UNIVERSAL_SAIGE_DIR GROUP_FILE.txt.gz CHROM [SPLICE_LABEL=spliceai2_gt040]
# The two AlphaMissense files of one chromosome, from brava_create_annot.py's group file and long table
# (GROUP_FILE.long.tsv.gz):
#   GROUP.flexrv_AM.txt.gz            the group file plus a 'score:AM' line per gene, for FlexRV (universal-saige step 2
#                                     --flexRVscore AM): pLoF 1.0, missense variants their AlphaMissense score
#   GROUP.am_moderate.saige_group.txt.gz   missense variants with AlphaMissense >= 0.906 (ACMG/AMP PP3 Moderate),
#                                     label am_moderate, for a separate SKAT-O run
# The scores come from universal-saige's 04_flexrv_groupfile.sh (the canonical FlexRV builder), on each entry's own
# transcript (--annoTable). Its AlphaMissense, MANE and GENCODE v39 files: bash download_resources.sh --alphamissense
# in universal-saige. Every missense-type label can hold a missense variant -- including the SpliceAI category (a
# missense variant with low REVEL/CADD and a splice call) -- so all three are passed to --missenseAnno, and the FlexRV
# mask in step 2 is then pLoF:damaging_missense:other_missense:SPLICE_LABEL.
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
( cd "$US" && bash 04_flexrv_groupfile.sh --group "$g" --chr "$c" --name AM \
    --missenseAnno damaging_missense,other_missense,$lab --nonMissense drop \
    --am "$am" --isoforms "$iso" --annoTable "$g.long.tsv.gz" --mane "$R/MANE.GRCh38.summary.txt.gz" \
    --gtf "$R/gencode.v39.annotation.gtf.gz" --provenance "$stem.flexrv_AM.provenance.tsv" --out "$stem.flexrv_AM.txt" )
gzip -f "$stem.flexrv_AM.txt"
python3 "$here/am_mask.py" "$g" "$stem.flexrv_AM.provenance.tsv" "$stem.am_moderate.saige_group.txt.gz"
gzip -f "$stem.flexrv_AM.provenance.tsv"
echo "wrote $stem.flexrv_AM.txt.gz and $stem.am_moderate.saige_group.txt.gz"
