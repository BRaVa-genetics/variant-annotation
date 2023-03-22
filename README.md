# Variant annotation for BRaVa
This repository applies a series of steps to annotate variants for gene-based testing in SAIGE.
1. Run VEP version 105 and the LOFTEE v1.04_GRCh38.
2. Run SpliceAI
3. Merge the resultant anntotations and use Hail to extract variant annotations according to recommendations. https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#
4. Annotate variants with MAC and MAF information in gnomAD

## VEP
- `.json` file used for VEP is here.

## SpliceAI
- This [fork](https://github.com/geertvandeweyer/SpliceAI) is reccomended as it enables batching of variant annotations on the GPU 
- CPU performance: ~1k predictions per hour
- GPU (A100 40GB) performance: ~700k predictions per hour

## Run SpliceAI (with Docker):
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
docker pull cmgantwerpen/spliceai_v1.3
docker run --gpus all cmgantwerpen/spliceai_v1.3:latest spliceai -I input.vcf -O output.vcf -R hg38.fa.gz -A grch38 -B 256 -T 32 
```
## Run SpliceAI (without Docker):
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

git clone https://github.com/geertvandeweyer/SpliceAI.git
cd SpliceAI
python setup.py install

pip install tensorflow

spliceai -I input.vcf -O output.vcf -R hg38.fa.gz -A grch38 -B 256 -T 32 
```
