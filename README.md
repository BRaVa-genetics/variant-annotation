# Variant annotation for BRaVa
This repository applies a series of steps to annotate variants for gene-based testing in SAIGE.
1. Run VEP version 105 and the LOFTEE v1.04_GRCh38.
2. Run SpliceAI
3. Merge the resultant anntotations and use Hail to extract variant annotations according to recommendations. https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#
4. Annotate variants with MAC and MAF information in gnomAD

## VEP
Install VEP version 105 via CLI as described [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) or summarized below:
```
# download ensembl repo
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep

# install vep version 105
git checkout release/105
perl INSTALL.pl
```


- `.json` file used for VEP is here.


## LOFTEE
Download LOFTEE version XX and dependencies via CLI:
```



```





## SpliceAI Annotations
- Using the fork https://github.com/geertvandeweyer/SpliceAI
- This fork enables batching of variant annotations on the GPU 
- CPU performance: ~1k predictions per hour
- GPU (A100 40GB) performance: ~700k predictions per hour
- Docker: https://hub.docker.com/r/cmgantwerpen/spliceai_v1.3/
