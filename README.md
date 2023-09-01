# BRaVa variant annotation
This repository contains all information and scripts required to generate annotation group-files ready to use for group tests in SAIGE-gene, and is split into three steps:

1. Run VEP version 105 with LOFTEE v1.04_GRCh38 (Docker/Singularity provided), and post-process the resultant VEP annotated vcf
2. Run SpliceAI
3. Run the Python BRaVa annotation script to extract variant annotations [according to recommendations](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#), and generate SAIGE annotation group-files

Note, before [**step 1**](#1-run-vep-version-105-with-loftee-v104_grch38), ensure that the VCF has split multiallelics, that variant IDs are of the form `CHROM:POS:REF:ALT`, and that you extract sites only vcf files from these, ready for VEP annotation to avoid huge I/O overheads:
```
bcftools view --drop-genotypes input.vcf.gz -O z -o sites_only_input.vcf.gz
```
where `input.vcf.gz` is your vcf file including genotype information that we wish to annotate, and `sites_only_output.vcf.gz` is the name of the compressed vcf file with all the genotype information from the samples removed.

If you don't have BCFtools (we'd be surprised though), go ahead and install it following the instructions [here](https://samtools.github.io/bcftools/howtos/install.html).

Importantly, we'll need to use a BCFtools plugin (split-vep) at the end of [**step 1**](#1-run-vep-version-105-with-loftee-v104_grch38). In order to use the BCFtools plugins, the environment variable `BCFTOOLS_PLUGIN` must be set and point to the correct location:

```
export BCFTOOLS_PLUGINS=/path/to/bcftools/plugins
```

(replacing `/path/to/bcftools/plugins` with the path to your BCFtools plugins folder). It may already be set within your compute environment, so make sure to check that first!

## 1. Run VEP version 105 with LOFTEE v1.04_GRCh38

Complete instructions and code are provided in our [VEP 105 LOFTEE repository](https://github.com/BRaVa-genetics/vep105_loftee). Briefly, the steps are:

- Download all of the required resources (cache etc)
- Run VEP in Docker/Singularity
- Post-process the file for input into [**step 3**](#3-run-the-python-brava-annotation-script-to-extract-variant-annotations)

Don't forget to run the post-processing step which removes variants with gnomAD [popmax](https://gnomad.broadinstitute.org/help/popmax) > 0.01 and uses the split-vep plugin for BCFtools (this simply splits multiple transcript annotations for a given variant across multiple lines, creating a table with a constant number of columns).

Following completion of [**step 1**](#1-run-vep-version-105-with-loftee-v104_grch38) for each of your vcfs, you should have a series of files containing the following columns:

| SNP_ID | GENE | LOF | MAX_AF | REVEL_SCORE | CADD_PHRED | CSQ | TRANSCRIPT | MANE_SELECT | CANONICAL | BIOTYPE |
| ------ | ---- | --- | ------ | ----------- | ---------- | --- | ---------- | ----------- | --------- | ------- | 

We will use this output, together with the output of [**step 2**](#2-run-spliceai) to define our annotations according to the [recommendations](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#). Now, let's move onto **step 2**:

## 2. Run SpliceAI
To run SpliceAI we provide two main options: a batched version of SpliceAI that uses GPUs for increased performance, and the original SpliceAI. For each of these options, **please follow our instructions in this repository** rather than those in the various READMEs within the repositories that we link to. This is to ensure that we align to the same gencode and ensembl versions as detailed in our [annotation recommendations](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#).

Note that the required [gencode.v39.ensembl.v105.annotation.txt.gz](https://github.com/BRaVa-genetics/variant-annotation/tree/main/data/SpliceAI) file is hosted in this repository.

### Batched version (recommended if GPUs are available!)
- This [fork](https://github.com/geertvandeweyer/SpliceAI) is highly recommended as it enables batching of variant annotations on the GPU 
- CPU performance: ~1k predictions per hour
- GPU (A100 40GB) performance: ~700k predictions per hour 🚀

#### SpliceAI Batching parameters
*When setting the batching parameters, be mindful of the system and GPU memory of the machine you are running the script on. Feel free to experiment, but some reasonable -T numbers would be 64/128/256. CPU memory is larger, and increasing -B might further improve performance.*

#### Running batched SpliceAI (with Docker):
First, navigate to the directory this README is located in. Then, run
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P ./data/
docker pull cmgantwerpen/spliceai_v1.3
docker run --gpus all cmgantwerpen/spliceai_v1.3:latest spliceai -I input.vcf -O output.vcf -R ./data/hg38.fa.gz -A ./data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz -B 4096 -T 256
```
where `input.vcf` is the name of your (uncompressed) sites only vcf file for splice annotation, and `output.vcf` is is your desired output filename of the resultant annotated uncompressed vcf.

#### Running batched SpliceAI (without Docker):
First, navigate to the directory this README is located in. Then, run
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P ./data/

git clone https://github.com/geertvandeweyer/SpliceAI.git
cd SpliceAI
python setup.py install

pip install tensorflow

spliceai -I input.vcf -O output.vcf -R ./data/hg38.fa.gz -A ./data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz -B 4096 -T 256
```
where `input.vcf` is the name of your (uncompressed) sites only vcf file for splice annotation and `output.vcf` is your desired output filename of the resultant annotated uncompressed vcf.

### Vanilla SpliceAI
#### Using pip
First, navigate to the directory this README is located in. Then, run
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P ./data/

pip install spliceai
pip install tensorflow

spliceai -I input.vcf -O output.vcf -R ./data/hg38.fa.gz -A ./data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz
```
where `input.vcf` is the name of your (uncompressed) sites only vcf file for splice annotation and `output.vcf` is your desired output filename of the resultant annotated uncompressed vcf.

#### Directly from Illumina's github repository
First, navigate to the directory this README is located in. Then, run
```
git clone https://github.com/Illumina/SpliceAI.git
cd SpliceAI
python setup.py install

pip install tensorflow

spliceai -I input.vcf -O output.vcf -R ./data/hg38.fa.gz -A ./data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz
```
where `input.vcf` is the name of your (uncompressed) sites only vcf file for splice annotation and `output.vcf` is your desired output filename of the resultant annotated uncompressed vcf.

Great, now you're here you should have a collection of SpliceAI annotated vcf files, and a collection of tab delimited text files munged from the output of VEP. We can now combine these and generate the annotations for input into gene-based testing in the next step:

## 3. Run the Python BRaVa annotation script to extract variant annotations

Pass your processed (tab-delimited) VEP file and your SpliceAI vcf file with the appropriate arguments to generate SAIGE group files ready for analysis in [universal-saige](https://github.com/BRaVa-genetics/universal-saige/).

If you used our instructions for [**step 1**](#1-run-vep-version-105-with-loftee-v104_grch38), then you can save a lot of typing and run:
```
python SAIGE_annotations/brava_create_annot.py -v vep_table -s spliceai_vcf -w output_file
```
where, `vep_table` is the path to the (tab-delimited) VEP annotations text file from the output of [**step 1**](#1-run-vep-version-105-with-loftee-v104_grch38) and `spliceai_vcf` is the path to annotated output file from spliceAI in **step 2**, as all of the defaults for the VEP columns are set at the names in the output of [**step 1**](#1-run-vep-version-105-with-loftee-v104_grch38).

Complete usage and requirements are in the [SAIGE_annotations folder](https://github.com/BRaVa-genetics/variant-annotation/tree/main/SAIGE_annotations). Note that this script allows for your own choice of column names for each of the required VEP columns.
