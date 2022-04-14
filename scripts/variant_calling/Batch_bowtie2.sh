#!/usr/bin/bash

# processes multiple fastq files in a directory

# Usage: ./Batch_bowtie2.sh <path/to/index-directory> <readsDirectory> <path/to/bowtie-shell-script.

indexpath=$1

readsdir=$2

pathtobowtie=$3

source /project/farman_uksr/miniconda3/etc/profile.d/conda.sh

conda activate py27

for f in `ls $readsdir/*_1*f*q*`; do pref=${f/_*/}; readID=${pref/*\//}; ./$pathtobowtie $indexpath $readsdir $readID; done

conda deactivate
