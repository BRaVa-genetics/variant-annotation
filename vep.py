#!/usr/bin/env python3

import hail as hl
import argparse
import os

from gnomad.utils.vep import process_consequences

def main(args):
    
    # input parameters
    vcf_path = args.vcf_path
    json_path = args.json_path
    log_path = args.log_path
    tmp_path = args.tmp_path
    out_prefix = args.out_prefix
    cpus = 1

    # init hail
    hl.init(
        log=log_path,
        default_reference="GRCh38",
        append=True,
        min_block_size=1,
        tmp_dir=tmp_path,
        local_tmpdir=tmo_path,
        master='local[%s]' % cpus # Required to prevent using more CPU resources than requested on HPC
    ) 

    # import vcf file
    mt = hl.import_vcf(vcf_path)
    
    # run vep and annotate variants
    mt = hl.vep(mt, json_path)

    # format variant to worst_csq_by_gene_canonical etc.
    mt = process_consequences(mt)
    
    # create hail table and write
    ht = mt.rows()
    ht.write(out_prefix + ".ht", overwrite=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_path', default=None, help='Path to input VCF')
    parser.add_argument('--json_path', default=None, help='Path to vep .json config file')
    parser.add_argument('--log_path', default=None, help='Path to log file')
    parser.add_argument('--tmp_path', default=None, help='Path to folder for temporary files')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for annotated output')
    args = parser.parse_args()

    main(args)



