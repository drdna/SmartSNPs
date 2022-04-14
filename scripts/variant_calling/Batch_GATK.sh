#!/usr/bin/bash

# processes multiple bowtie alignment directories

# Usage: ./Batch_GATK.sh <path/to/index-directory> <readsDirectory> <path/to/GATK-shell-script.

indexpath=$1

readsdir=$2

pathtoGATK=$3

source /project/farman_uksr/miniconda3/etc/profile.d/conda.sh

conda activate py27

for f in `ls $readsdir/*_1*f*q*`; do pref=${f/_*/}; readID=${pref/*\//}; ./$pathtoGATK $indexpath $readID; done

conda deactivate
