# BRaVa SAIGE annotation group files

*Step 4 of the [BRaVa annotation pipeline](https://github.com/BRaVa-genetics/variant-annotation).*

`scripts/brava_create_annot.py` merges the VEP table (step 1) with SpliceAI2 scores (step 2; or SpliceAI v1) and
writes a SAIGE group file with BRaVa's labels, plus a long table with one row per (gene, variant) and every field its
label came from. The rules are in the [top-level README](../README.md#the-labels) and the script's docstring.

## Requirements

Python 3, standard library only. `python3 scripts/brava_create_annot.py --selftest` runs 31 checks on built-in
examples (the rules at each cut-off, REVEL list formats, strand and gene-body handling of SpliceAI2, the MANE Select
row winning over a canonical one, and the refusals below).

## Usage

```bash
python3 brava_create_annot.py -v VEP_TABLE --spliceai2 SCORES.tsv.gz [SCORES.tsv.gz ...] --gtf gencode.v39.annotation.gtf.gz -w OUT.txt.gz
python3 brava_create_annot.py -v VEP_TABLE -s SPLICEAI_V1.vcf -w OUT.txt.gz
```

## Arguments

* `-v/--vep`: the VEP table from step 1 (whitespace-delimited, header line).
* `--spliceai2`: SpliceAI2 score tables (columns chrom, pos, ref, alt, strand, ..., score last; a header line is
  skipped): the precomputed lookups and the GPU-scored chunks of step 2. Needs `--gtf`.
* `--gtf`: GENCODE v39 GTF. Each entry takes the SpliceAI2 score on its gene's strand, and only when the variant lies
  inside the gene (otherwise the score is kept in `DS_MAX` with `SPLICEAI_SOURCE` = `spliceai2_outside_gene` and
  gives no splice route).
* `-s/--spliceai`: a SpliceAI v1 VCF instead of SpliceAI2 (label `spliceai_gt05`). `chr1` and `1` both work.
* `-w/--out_file`: the group file (gzipped if it ends in `.gz`); the long table goes to `OUT_FILE.long.tsv.gz`.
* `--cadd_indels`: CADD output for indels (step 3), used where VEP has no CADD.
* `--allow_missing_loftee`: write the files even when LoF-type entries have no LOFTEE call (see below).
* `--vep_snp_id_col`, `--vep_gene_col`, `--vep_lof_col`, `--vep_revel_col`, `--vep_cadd_phred_col`,
  `--vep_consequence_col`, `--vep_transcript_col`, `--vep_mane_select_col`, `--vep_canonical_col`,
  `--vep_biotype_col`: column names in the VEP table, if yours differ from step 1's.

## It refuses to write when

* **a stop_gained / frameshift / splice donor or acceptor entry has no LOFTEE call.** LOFTEE calls every such entry HC
  or LC, so a missing call means LOFTEE failed on that variant (VEP's `--fork` with LOFTEE's conservation database
  did this at random; [vep105_loftee](https://github.com/BRaVa-genetics/vep105_loftee#why-no-conservation-database)).
  Re-run VEP on those variants and merge;
* **no SpliceAI score matches any entry**: almost always a mismatch of variant IDs, chromosome names or genome build.

## Output

* the group file: per gene a `GENE var ...` line and a `GENE anno ...` line, variants in the order of the VEP table;
* the long table: `ID GENE LOF REVEL_SCORE CADD_PHRED CSQ TRANSCRIPT MANE_SELECT CANONICAL BIOTYPE DS_MAX
  SPLICEAI_SOURCE ANNOTATION`. It is the `--annoTable` input of universal-saige's FlexRV builder
  (`scripts/alphamissense_files.sh`).

## Checked against the October 2026 release

On chromosomes 21, 22 and X, the script reproduces the label, LoF call, consequence, transcript and SpliceAI2 score of
all 661,493 entries of the BRaVa October 2026 group files that were annotated from scratch (re-annotated variants of
six biobanks), with no entry missing.
