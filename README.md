# BRaVa Variant Annotation
This repository contains all information and scripts required to generate SAIGE annotations for group tests.

1. Run VEP version 105 w/ LOFTEE v1.04_GRCh38 (Docker/Singularity provided)
2. Run SpliceAI (batched version)
3. Run the Python BRaVa annotation script to extract variant annotations [according to recommendations](https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#) and generate SAIGE group files.

## 1. Run VEP version 105 w/ LOFTEE v1.04_GRCh38

Complete instructions and code provided in the [vep105_loftee repo](https://github.com/BRaVa-genetics/vep105_loftee)

- Download all required resources (cache etc...)
- Run VEP in Docker/Singularity
- Generate the "worst consequence by gene (MANE SELECT canonical)" variants from the VEP output

## 2. Run SpliceAI (batched version)
- This [fork](https://github.com/geertvandeweyer/SpliceAI) is highly recommended as it enables batching of variant annotations on the GPU 
- CPU performance: ~1k predictions per hour
- GPU (A100 40GB) performance: ~700k predictions per hour 🚀
- [gencode.v39.ensembl.v105.annotation.txt.gz is hosted in this repo](https://github.com/BRaVa-genetics/variant-annotation/tree/main/data/SpliceAI)

### SpliceAI Batching parameters
*When setting the batching parameters, be mindful of the system and gpu memory of the machine you are running the script on. Feel free to experiment, but some reasonable -T numbers would be 64/128. CPU memory is larger, and increasing -B might further improve performance.*

### Running SpliceAI (with Docker):
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P ./data/
docker pull cmgantwerpen/spliceai_v1.3
docker run --gpus all cmgantwerpen/spliceai_v1.3:latest spliceai -I input.vcf -O output.vcf -R ./data/hg38.fa.gz -A ./data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz -B 4096 -T 256
```
### Running SpliceAI (without Docker):
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P ./data/

git clone https://github.com/geertvandeweyer/SpliceAI.git
cd SpliceAI
python setup.py install

pip install tensorflow

spliceai -I input.vcf -O output.vcf -R ./data/hg38.fa.gz -A ./data/SpliceAI/gencode.v39.ensembl.v105.annotation.txt.gz -B 4096 -T 256
```

## 2. Run the Python BRaVa annotation script to extract variant annotations

Complete instructions and code provided in the [SAIGE-annotations-for-BRaVa repo](https://github.com/BRaVa-genetics/variant-annotation-python/tree/main)

- Pass your processed VEP file and your SpliceAI vcf file with the appropriate arguments
- Generate final SAIGE group files ready for analysis in [universal-saige](https://github.com/BRaVa-genetics/universal-saige/) ✅
