# BRaVa SAIGE Annotation

*This repository is [part 3](https://github.com/BRaVa-genetics/variant-annotation#3-run-the-python-brava-annotation-script-to-extract-variant-annotations) of the full BRaVa annotation pipeline. Instructions for the full pipeline are [here](https://github.com/BRaVa-genetics/variant-annotation).*

The python script in this repository merges VEP and SpliceAI information to generate a SAIGE annotation groupfile with annotations defined according to BRaVa's [guidelines](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#). The required inputs for the annotations by VEP and SpliceAI are table format files (tab-delimited).

The bundled shell script (`scripts/brava_create_annot.sh`) provides an example of running these commands on an HPC cluster as an array job.

## Requirements

This script requires `numpy` and `pandas` to be installed in your Python environment.

## Usage

```bash
python brava_create_annot.py -c chr -v vep_table -s spliceai_vcf -w output_file --vep_snp_id_col snp_id_col --vep_gene_col gene_id_col --vep_lof_col lof_col --vep_max_pop_col gnomad_popmax_AF_col --vep_revel_col revel_col --vep_cadd_phred_col cadd_col --vep_consequence_col csq_col --vep_canonical_col canonical_col --vep_biotype_col biotype_col --vep_mane_select_col mane_select_col
```

Note that if you followed our [**step 1**](https://github.com/BRaVa-genetics/variant-annotation#1-run-vep-version-105-with-loftee-v104_grch38) commands, you won't need to use any arguments other than `-c`, `-v`, `-s` and `-w`, as these are set at the defaults.

## Arguments

* `-c/--chr`: chromosome to process.
* `-v/--vep`: file path to the (tab-delimited) VEP annotated table ([**step 1**](https://github.com/BRaVa-genetics/variant-annotation#1-run-vep-version-105-with-loftee-v104_grch38) output)
* `-s/--spliceai`: file path to the SpliceAI annotated vcf file ([**step 2**](https://github.com/BRaVa-genetics/variant-annotation#2-run-spliceai) output).
* `-w/--out_file`: output file where the results will be saved
* `--vep_snp_id_col`: SNP ID (chr:pos:ref:alt) column in the VEP annotated table
* `--vep_gene_col`: Gene ID column in the VEP annotated table
* `--vep_lof_col`: LoF column in the VEP annotated table
* `--vep_max_pop_col`: gnomAD_maxAF column in the VEP annotated table
* `--vep_revel_col`: REVEL column in the VEP annotated table
* `--vep_cadd_phred_col`: Phred scaled CADD score column in the VEP annotated table
* `--vep_consequence_col`: consequence column in the VEP annotated table
* `--vep_canonical_col`: canonical transcript column in the VEP annotated table
* `--vep_biotype_col`: biotype column in the VEP annotated table
* `--vep_mane_select_col`: MANE select column in the VEP annotated table (not that this the column that annotates the name of the transcript, not a boolean identifier of whether the annotated transcript is MANE select (annotated using the `--mane-select` flag))

