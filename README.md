# Variant annotation for BRaVa
This repository applies a series of steps to annotate variants for gene-based testing in SAIGE.
1. Run VEP version 105 and the LOFTEE v1.04_GRCh38.
2. Run SpliceAI
3. Merge the resultant anntotations and use Hail to extract variant annotations according to recommendations. https://docs.google.com/document/d/11Nnb_nUjHnqKCkIB3SQAbR6fl66ICdeA-x_HyGWsBXM/edit#
4. Annotate variants with MAC and MAF information in gnomAD


## Installing Hail-compatible VEP with conda
- Install the conda env below directly from the .yml file:
```
conda env create -f bravavep.yml

```

- If this does not work, then code for creating the .yml and building baseline VEP v105 can be found here:
```
# ensure we have the right channels
conda config --append channels conda-forge 
conda config --append channels bioconda
# create environemnt 
conda create -n brava-vep-v3
conda install -c bioconda ensembl-vep=105.0
pip install hail notebook
# install zlib (only conda-forge has v2.202-pl5321h166bdaf_0)
conda install -c conda-forge perl-compress-raw-zlib 
```

- Download and install human vep cache:
```
my_vep_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/brava/vep_test"
vep_install -a cf -s homo_sapiens -y GRCh38 -c ${my_vep_dir} --CONVERT
(ls ${my_vep_dir}/homo_sapiens/105_GRCh38  && echo "SUCCESS") || echo "FAILED"
```

- Test if VEP can be run from the command line by typing 'vep'.
- Test if VEP can be run through hail by opening up a python session with 'python3' and entering:
```
import hail as hl
from gnomad.utils.vep import process_consequences

hl.init(
        log="test.log",
        default_reference="GRCh38",
        append=True,
        min_block_size=1,
        tmp_dir='/well/lindgren-ukbb/projects/ukbb-11867/flassen/spark',
        local_tmpdir='/well/lindgren-ukbb/projects/ukbb-11867/flassen/spark',
        master='local[%s]' % 1 # Required to prevent using more CPU resources than requested
    )

mt = hl.read_matrix_table("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/unphased/wes/prefilter/200k/ukb_split_wes_200k_chr21_parents.mt")
vep_mt = hl.vep(mt, "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/configs/vep_test.json")
>>> vep_mt.vep.most_severe_consequence
<StringExpression of type str>
>>> vep_mt.vep.most_severe_consequence.show()
+----------------+------------+--------------------------------------+
| locus          | alleles    | <expr>                               |
+----------------+------------+--------------------------------------+
| locus<GRCh38>  | array<str> | str                                  |
+----------------+------------+--------------------------------------+
| chr21:10413617 | ["C","T"]  | "non_coding_transcript_exon_variant" |
| chr21:10413618 | ["G","A"]  | "non_coding_transcript_exon_variant" |
| chr21:10413624 | ["C","T"]  | "non_coding_transcript_exon_variant" |
| chr21:10413627 | ["C","G"]  | "non_coding_transcript_exon_variant" |
| chr21:10413629 | ["C","T"]  | "non_coding_transcript_exon_variant" |
| chr21:10413631 | ["C","T"]  | "non_coding_transcript_exon_variant" |
| chr21:10413633 | ["C","T"]  | "non_coding_transcript_exon_variant" |
| chr21:10413634 | ["T","G"]  | "non_coding_transcript_exon_variant" |
| chr21:10413636 | ["G","T"]  | "non_coding_transcript_exon_variant" |
| chr21:10413638 | ["T","A"]  | "non_coding_transcript_exon_variant" |
+----------------+------------+--------------------------------------+
```





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
