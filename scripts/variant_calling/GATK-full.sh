#!/bin/bash

# GATK version currently in use: 4.1.4.1

# Usage GATK-full.sh <path/to/index-directory> <readsID (prefix only)>


# specify path to conda script

source /project/farman_uksr/miniconda3/etc/profile.d/conda.sh

# specify index information

indexpath=$1

indexdir=${indexpath/*\//}

indexname=${indexdir/_*/}

fasta=$indexpath/${indexname}".fasta"


# give reads prefix (base ID)

readsprefix=$2


# identify bowtie output directory

indir=${indexname}"_"${readsprefix}"_ALIGN"


# index reference fasta file

conda activate py27

samtools faidx $fasta

conda deactivate


# activate gatk environment 

conda activate gatk


# Create dict and fai for reference

# remove existing dictionary

rm ${fasta/fasta/dict}

#create dictionary

java -jar /project/farman_uksr/miniconda3/share/picard-2.21.6-0/picard.jar CreateSequenceDictionary R=$fasta O=${fasta/fasta/dict}

# Call haplotype

/project/farman_uksr/gatk-4.1.4.1/gatk --java-options "-Xmx35g" HaplotypeCaller \
        --native-pair-hmm-threads 16 \
	-R $fasta \
        -ploidy 1 \
        -I $indir/accepted_hits_sortedRG.bam \
        --emit-ref-confidence GVCF \
        -O $indir/${readsprefix}.vcf \

# Determine genotypes

/project/farman_uksr/gatk-4.1.4.1/gatk --java-options "-Xmx35g" GenotypeGVCFs \
 -R $fasta \
 -V $indir/$readsprefix.vcf \
 --output $indir/${readsprefix}_genotyped.vcf


# Filter SNPs only

/project/farman_uksr/gatk-4.1.4.1/gatk SelectVariants \
 -R $fasta \
 -select-type SNP \
 -V $indir/${readsprefix}_genotyped.vcf \
 --output $indir/${readsprefix}_genotyped-snps.vcf


vcftools --vcf $indir/${readsprefix}_genotyped.snp-only.filters.vcf \
        --remove-filtered-all --recode \
        --out $indir/${readsprefix}_genotyped.snp-only.filters.subset.vcf

conda deactivate
