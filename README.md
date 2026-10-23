# BRaVa variant annotation

Everything needed to build BRaVa's SAIGE annotation group files from your own variant sites, ready for gene-based
testing with [universal-saige](https://github.com/BRaVa-genetics/universal-saige).

> [!IMPORTANT]
> **What's new (October 2026)**
> 1. **SpliceAI has its own category.** SpliceAI no longer routes a variant into `damaging_missense`
>    (previously DS >= 0.2). A splice prediction at the ACMG/AMP PP3 Moderate level -- the same evidence level as
>    REVEL >= 0.773 / CADD >= 28.1 -- is labelled `spliceai2_gt040` (SpliceAI2 >= 0.40, the default) or
>    `spliceai_gt05` (SpliceAI v1 >= 0.5).
> 2. **SpliceAI2** (Illumina, 2026) is the default splice model; SpliceAI v1 stays available. One run uses one model.
> 3. **Short label names**: `pLoF`, `damaging_missense`, `spliceai2_gt040`, `other_missense`, `synonymous`,
>    `non_coding` (were `damaging_missense_or_protein_altering`, `other_missense_or_protein_altering`).
> 4. **dbNSFP 5.4a** (REVEL, CADD 1.7) for everyone, and **no LOFTEE conservation database**: with VEP's `--fork` it
>    made LOFTEE drop LoF calls at random (see [vep105_loftee](https://github.com/BRaVa-genetics/vep105_loftee#why-no-conservation-database)).
>    The annotation script now refuses to write group files when a LoF-type variant has no LOFTEE call.
> 5. **AlphaMissense files** for extra runs: an `am_moderate` mask (AlphaMissense >= 0.906) and a FlexRV file on the
>    FlexRV paper's scale (`flexrv_AMPCT`: each missense variant's within-gene AlphaMissense percentile).
> 6. `brava_create_annot.py` needs only Python 3 (no numpy / pandas) and has a `--selftest`.

### Contents
* [The labels](#the-labels)
* [Before you start](#before-you-start)
* [1. Run VEP 105 with LOFTEE](#1-run-vep-105-with-loftee)
* [2. SpliceAI2 scores (default)](#2-spliceai2-scores-default)
* [2b. SpliceAI v1 instead](#2b-spliceai-v1-instead)
* [3. CADD for indels (optional)](#3-cadd-for-indels-optional)
* [4. Build the group files](#4-build-the-group-files)
* [5. AlphaMissense files](#5-alphamissense-files)
* [6. Run SAIGE](#6-run-saige)

## The labels

One label per (gene, variant), on the gene's MANE Select transcript, else its canonical one (protein-coding
transcripts only). The first rule that applies wins:

| label | rule |
|---|---|
| `pLoF` | LOFTEE high confidence |
| `damaging_missense` | missense-type and (REVEL >= 0.773 or CADD >= 28.1); or LOFTEE low confidence |
| `spliceai2_gt040` | SpliceAI2 >= 0.40 (with SpliceAI v1: `spliceai_gt05`, DS >= 0.5) |
| `other_missense` | other missense-type variants |
| `synonymous` | synonymous with SpliceAI2 < 0.10 (v1: < 0.2) |
| `non_coding` | everything else, including synonymous with SpliceAI2 0.10-0.40 (v1: 0.2-0.5) |

Missense-type: `missense_variant`, `inframe_insertion`, `inframe_deletion`, `start_lost`, `stop_lost`,
`protein_altering_variant`, `transcript_amplification`. A consequence belongs to a class when any of its
`&`-joined VEP terms does. Variants with a gnomAD v2 exome popmax > 0.01 are removed in step 1. A variant with no
SpliceAI score has no splice route. Why 0.40 and 0.5: both are the PP3 Moderate cut-offs from an ACMG/AMP
calibration of each model against ClinVar (Illumina's README pairs SpliceAI2 0.25 with v1's 0.5, but 0.25 sits at
Supporting).

## Before you start

Split multiallelic variants, set every variant ID to `CHROM:POS:REF:ALT`, and keep a sites-only VCF per chromosome:
the commands are in [vep105_loftee](https://github.com/BRaVa-genetics/vep105_loftee#pre-processing). Needs `bcftools`
and `tabix`/`bgzip`; step 2 needs a CUDA GPU for the variants SpliceAI2's precomputed scores do not cover.

## 1. Run VEP 105 with LOFTEE

Follow [vep105_loftee](https://github.com/BRaVa-genetics/vep105_loftee): download the resources (VEP cache, LOFTEE
files, dbNSFP 5.4a), run VEP in Docker or Singularity, **check LOFTEE's warnings count is 0**, then post-process
(gnomAD popmax filter and `bcftools +split-vep`). You end with one table per chromosome with the columns

| SNP_ID | GENE | LOF | REVEL_SCORE | CADD_PHRED | CSQ | TRANSCRIPT | MANE_SELECT | CANONICAL | BIOTYPE |
| ------ | ---- | --- | ----------- | ---------- | --- | ---------- | ----------- | --------- | ------- |

## 2. SpliceAI2 scores (default)

SpliceAI2 (github.com/Illumina/SpliceAI2) is licensed by Illumina for non-commercial use by non-commercial
organisations only; its model weights and precomputed scores are gated on Hugging Face (accept the terms, then put a
read token in `~/.cache/huggingface/token`). Any publication using it must cite Illumina and Jaganathan K, Chen J,
Liu X, Zhang Y, et al., *A unified framework for quantitative splicing and transcript prediction* (bioRxiv 2026,
DOI 10.64898/2026.10.08.757673). Do not share the scores themselves outside your analysis.

SpliceAI2 scores each position once per strand, not per gene: an entry takes the score on its gene's GENCODE v39
strand, and only when the variant lies inside the gene (step 4 does this). Most variants are already scored in
Illumina's precomputed release (all SNVs and ~150M population indels in protein-coding genes); the rest are scored
on a GPU. Scripts are in [SpliceAI2/scripts](SpliceAI2/scripts); the GENCODE v39 GTF is
`resources/gencode.v39.annotation.gtf.gz` in universal-saige (`bash download_resources.sh --alphamissense`).
```
# keys: the (variant, gene strand) pairs your entries need
python3 SpliceAI2/scripts/make_keys.py gencode.v39.annotation.gtf.gz spliceai2/keys out/*_processed.txt
# precomputed scores: one task per chromosome and kind (24 SNV, 23 indel; there is no chrY indel file);
# streams the 348 GB release without storing it, and lists the keys it lacks
bash SpliceAI2/scripts/lookup_chr.sh chr21 snv spliceai2/keys spliceai2/scores
bash SpliceAI2/scripts/lookup_chr.sh chr21 indel spliceai2/keys spliceai2/scores
# the missing keys worth scoring (inside a transcript SpliceAI v1 scores, within 250 bp of an exon), in chunks
python3 SpliceAI2/scripts/gpu_set.py spliceai2/scores gencode.v39.annotation.gtf.gz \
    data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz spliceai2/gpu
# once: SpliceAI2, its models and hg38 in a work directory (~7 GB)
bash SpliceAI2/scripts/setup.sh sai2
# controls: 300 variants the release already scores, re-scored with every chunk; a chunk is kept only if they agree
python3 SpliceAI2/scripts/control.py pick spliceai2/scores spliceai2/gpu/controls.tsv
# on a GPU node, per chunk (A100, NPROC=4: 1,278 variants in 7 min 20 s including model loading)
bash SpliceAI2/scripts/score.sh sai2 spliceai2/gpu/chunk_000.tsv spliceai2/gpu/controls.tsv
```

## 2b. SpliceAI v1 instead

To use SpliceAI v1 (DS >= 0.5, label `spliceai_gt05`), score your sites with the GENCODE v39 / Ensembl 105
annotation in this repository ([data/SpliceAI](data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz)) and UCSC
hg38. The batched fork (github.com/geertvandeweyer/SpliceAI) runs ~700k predictions/hour on an A100:
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P ./data/
docker pull cmgantwerpen/spliceai_v1.3
docker run --gpus all cmgantwerpen/spliceai_v1.3:latest spliceai -I input.vcf -O output.vcf -R ./data/hg38.fa.gz \
    -A ./data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz -B 4096 -T 256
```
or, without Docker, `pip install spliceai tensorflow` and the same `spliceai` command (without `-B` / `-T`), ~1k
predictions/hour on a CPU. `input.vcf` is your uncompressed sites-only VCF.

## 3. CADD for indels (optional)

VEP's dbNSFP plugin gives CADD for SNVs only. To score in-frame indels as well, run CADD on them and pass its output
(columns Chrom, Pos, Ref, Alt, RawScore, PHRED) to step 4 with `--cadd_indels`; it is used only where VEP has no
CADD. Use CADD v1.7 on GRCh38, the version dbNSFP 5.4a carries for SNVs, following the
[CADD-scripts](https://github.com/kircherlab/CADD-scripts) README (its GRCh38 annotations are ~200 GB). CADD wants
chromosome names without `chr`:
```
bcftools view -v indels file.vcf.gz -Ov | sed 's/^chr//' > file_indels.vcf
```

## 4. Build the group files

```
python3 SAIGE_annotations/scripts/brava_create_annot.py --selftest      # once: 31 checks, no data needed
python3 SAIGE_annotations/scripts/brava_create_annot.py \
    -v out/sites_only_output_chr21_vep.gnomad_popmax_0.01_processed.txt \
    --spliceai2 spliceai2/scores/chr21.snv.tsv.gz spliceai2/scores/chr21.indel.tsv.gz spliceai2/gpu/chunk_*.scores.tsv.gz \
    --gtf gencode.v39.annotation.gtf.gz \
    -w brava.chr21.saige_group.txt.gz
```
(with SpliceAI v1: `-s output.vcf` in place of `--spliceai2 ... --gtf ...`). It writes the group file and
`brava.chr21.saige_group.txt.gz.long.tsv.gz`, one row per (gene, variant) with every field its label came from.
Details and all options: [SAIGE_annotations](SAIGE_annotations).

## 5. AlphaMissense files

Per chromosome, from the group file and its long table, using universal-saige's FlexRV builder:
```
bash SAIGE_annotations/scripts/alphamissense_files.sh /path/to/universal-saige brava.chr21.saige_group.txt.gz 21
```
writes two files:

* `brava.chr21.saige_group.flexrv_AMPCT.txt.gz`: FlexRV weights on the FlexRV paper's scale (Schwartzentruber et al.,
  bioRxiv 2025.11.04.686320). Each score is computed as follows:

  | entry | score |
  |---|---|
  | pLoF | 1.0 |
  | missense | percentile of its raw AlphaMissense score among all possible missense variants of the gene's transcript: `(#below + (#equal+1)/2) / N` |
  | protein-altering non-missense (in-frame indel, LOFTEE-LC LoF, start/stop lost) | the gene's mean, after ranking |
  | no protein change | dropped |
  | missense in a paralog copy AlphaMissense does not score (identical protein) | its protein-level score for the same amino-acid change, then as above |
  | missense in a gene AlphaMissense does not score | 0.5 (the expected within-gene percentile) |

  The transcript is MANE Select, else canonical, in the version GENCODE v39 holds. The full mapping is in
  universal-saige's README, step 4.
* `brava.chr21.saige_group.am_moderate.saige_group.txt.gz`: missense variants with raw AlphaMissense >= 0.906 (the
  PP3 Moderate cut-off), each scored from the same AlphaMissense row as its `flexrv_AMPCT` entry. Both files come from
  one run.

With SpliceAI v1, add `spliceai_gt05` as a fourth argument. `KEEP_RAW=1` also keeps `flexrv_AM.txt.gz` (the raw-score
FlexRV file of the earlier release).

A fourth run uses the scallion weights. Its `flexrv_SCALLION.txt.gz` files are built centrally from the scallion
authors' within-gene percentiles, with the same variants and rules as `flexrv_AMPCT`, and come with the release.

## 6. Run SAIGE

With [universal-saige](https://github.com/BRaVa-genetics/universal-saige) step 2, four runs per chromosome and
phenotype (its [step 2 template](https://github.com/BRaVa-genetics/universal-saige/blob/main/templates/step_2_template.sh)
has them written out):

| run | file | `--annotations` |
|---|---|---|
| SAIGE-GENE+ (SKAT-O) | `*.saige_group.txt.gz` | `pLoF,damaging_missense,pLoF:damaging_missense,pLoF:damaging_missense:spliceai2_gt040,synonymous,other_missense,pLoF:damaging_missense:spliceai2_gt040:other_missense:synonymous` |
| SKAT-O, AlphaMissense >= 0.906 | `*.am_moderate.saige_group.txt.gz` | `am_moderate`, with `--singleInGroup FALSE` |
| FlexRV, AlphaMissense percentile weights | `*.flexrv_AMPCT.txt.gz`, with `--flexRVscore AMPCT` | `pLoF:damaging_missense:other_missense:spliceai2_gt040` |
| FlexRV, scallion weights | `*.flexrv_SCALLION.txt.gz`, with `--flexRVscore SCALLION` | `pLoF:damaging_missense:other_missense:spliceai2_gt040` |

With SpliceAI v1, use `spliceai_gt05` in place of `spliceai2_gt040`.
