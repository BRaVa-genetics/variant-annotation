# Variant annotation for BRaVa
This repository applies a series of steps to annotate variants for gene-based testing in SAIGE.
1. Run VEP version 105 and the LOFTEE v1.04_GRCh38.
2. Run SpliceAI
3. Merge the resultant anntotations and use Hail to extract variant annotations according to recommendations. https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#
4. Annotate variants with MAC and MAF information in gnomAD


## VEP (With docker)

## VEP (without docker)


### Installation
Install VEP version 105 via CLI as described [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) or summarized below:
```
# download ensembl repo
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep

# pull version 105
git checkout release/105
perl INSTALL.pl
```

- `.json` file used for VEP is here.

###  dbNSFP plguin installation
Guide adopted from [here](https://sites.google.com/site/jpopgen/dbNSFP).


### LOFTEE plugin installation
- Download LOFTEE:
```
git clone https://github.com/konradjk/loftee.git

```
- Download `human_ancestor_fa` dependency:
```
curl -O https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz
curl -O https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai
curl -O https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.gzi
```

- Download `conservation_file` dependency:
```
curl -O https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz
```








## SpliceAI Annotations
- Using the fork https://github.com/geertvandeweyer/SpliceAI
- This fork enables batching of variant annotations on the GPU 
- CPU performance: ~1k predictions per hour
- GPU (A100 40GB) performance: ~700k predictions per hour
- Docker: https://hub.docker.com/r/cmgantwerpen/spliceai_v1.3/
