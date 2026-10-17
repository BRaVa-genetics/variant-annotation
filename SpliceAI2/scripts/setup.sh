#!/bin/bash
# One-off: SpliceAI2 (github.com/Illumina/SpliceAI2; PyTorch, a CUDA GPU is required) in a conda prefix, plus the two
# model folders it averages over (Hugging Face illumina-ai/SpliceAI2, gated: accept its terms, then put a read token
# in ~/.cache/huggingface/token). Versions pinned to those the models were trained with (model README: Python 3.13.2,
# PyTorch 2.6.0+cu124, Lightning 2.6.1, NumPy 2.4.3, pandas 3.0.1, pyfaidx 0.9.0.3). ~7 GB.
# usage: setup.sh WORK_DIR      (then: WORK_DIR/sai2env, WORK_DIR/models/SpliceAI2, WORK_DIR/ref/hg38.fa)
set -euo pipefail
W=$(mkdir -p "$1" && cd "$1" && pwd); cd "$W"
export PIP_CACHE_DIR=$W/.cache/pip CONDA_PKGS_DIRS=$W/.cache/conda HF_HOME=$W/.cache/hf
mkdir -p "$PIP_CACHE_DIR" "$CONDA_PKGS_DIRS" "$HF_HOME" models ref
cp -n ~/.cache/huggingface/token "$HF_HOME/token" 2>/dev/null || true
ENV=$W/sai2env
if [[ ! -x $ENV/bin/python ]]; then
  CONDA=conda; command -v mamba >/dev/null && CONDA=mamba
  $CONDA create -y -p "$ENV" -c conda-forge python=3.13.2 pip
fi
"$ENV/bin/pip" install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
"$ENV/bin/pip" install lightning==2.6.1 numpy==2.4.3 pandas==3.0.1 pyfaidx==0.9.0.3 "huggingface_hub[cli]" spliceai2
"$ENV/bin/hf" download illumina-ai/SpliceAI2 --local-dir models/SpliceAI2
ls -l models/SpliceAI2/spliceai2_multispecies_13M_0{1,2}/model.ckpt
if [[ ! -s ref/hg38.fa ]]; then
  curl -sSfL http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gzip -dc > ref/hg38.fa
fi
"$ENV/bin/python" -c "import pyfaidx; pyfaidx.Fasta('ref/hg38.fa'); import torch, spliceai2; print('torch', torch.__version__, 'cuda build', torch.version.cuda)"
echo "setup done in $W"
