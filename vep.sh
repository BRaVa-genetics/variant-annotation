#!/usr/bin/env bash

conda activate brava-vep-v3

readonly tmp_dir="tmp/spark"
readonly hail_script="vep.py"
readonly in_json="vep.json"
readonly in_vcf="test/example.vcf"
readonly out_prefix="test/annotated_example"

mkdir -p ${tmp_dir}

python3 ${hail_script} \
     --vcf_path "${in_vcf}" \
     --json_path "${in_json}" \
     --tmp_dir "${tmp_dir}" \
     --out_prefix "${out_prefix}"





