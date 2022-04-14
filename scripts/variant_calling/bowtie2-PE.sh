#!/bin/bash

# Bowtie 2 version currently in use: 2.2.5

# Usage: bowtie2 <index-dir> <reads-dir> <reads-prefix>

source /project/farman_uksr/miniconda3/etc/profile.d/conda.sh

conda activate py27


# path to index directory

indexpath=$1


# path to directory containing reads

readsdir=$2


# specify reads prefix (e.g. genome1.fastq.gz prefix = genome1)

reads=$3	


# create output directory with format: indexID_readsID_ALIGN

indexdir=${indexpath/*\//}

indexname=${indexdir/_*/}

echo indexdir = $indexdir

echo indexname = $indexname

outdir=${indexname}_${reads}_ALIGN

mkdir $outdir


# run alignment

bowtie2 --threads 16 --very-sensitive-local --phred33 --no-unal -x $indexpath/$indexname \
 -1 $readsdir/${reads}_*1.f*q* -2 $readsdir/${reads}_*2.f*q* | sed 's/#0\/[1-2]//' | samtools view -bS - > $outdir/accepted_hits.bam


# sort output bamfile and remove unsorted data

samtools sort $outdir/accepted_hits.bam -o $outdir/accepted_hits_sorted.bam

rm $outdir/accepted_hits.bam


#conda deactivate


# add read group info to bam header

# define read group info

o=${reads}


# activate gatk environment

conda activate gatk


# picard version:2.21.6-0

java -jar /project/farman_uksr/miniconda3/share/picard-2.21.6-0/picard.jar AddOrReplaceReadGroups I=$outdir/accepted_hits_sorted.bam \
O=$outdir/accepted_hits_sortedRG.bam RGID=$o RGSM=$o RGLB=$o RGPI=50 RGPL=illumina RGPU=unit1

conda deactivate


# remove data without read group header info

rm $outdir/accepted_hits_sorted.bam


# index the bamfile

conda activate py27

# samtools version currently running: 1.9

samtools index $outdir/accepted_hits_sortedRG.bam

conda deactivate

# done
