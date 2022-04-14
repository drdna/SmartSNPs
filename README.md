# SmartSNPs

SmartSNPs is a postprocessor that uses hard and fast rules to filter SNPs that have been called in haploid fungal genomes. Variant calls generated based on read alignments to a reference genome are typically contaminated with large numbers of "illegal" calls. These occur because reads from repeated regions of the genome sometimes fail to align to genomic repeats in the reference genome, owing to significant sequence divergence caused by the gradual accumulation of mutations over time or, in fungi, by the action of Repeat-Induced Point mutation (RIP) - a mutagenec process that causes repeated sequences to esxperience massive mutational sweeps in a single sexual cycle. As a result, variants are often called in repeated regions of the genome that are not recognized as being repeated due to the unexpectedly low read coverage across each copy. Note: with highly diverged repeats, read alignment is stochastic because it depends not only on sequence divergence but also the positions of mismatches in the sequence reads.

Phylogenetic studies absolutely require that all each position in an alignment should be from a single allele. Therefore, unless one can phase the variants (a challenging proposition at best), SNP calls in repeated sequences absolutely should not be used for phylogenetic inference. In fungi, a second reason to filter SNPs in repeated sequences is because these mutations may be the direct or indirect products of RIP and, therefore, disobey fundamental evolutionary assumptions. Here it is important to recognize that if a SNP is found to occur in a sequence that is repeated in one genome, one can no longer be certain that the SNP calls in other genomes are truly allelic with the reference. For this reason, it is safest to purge calls that are repeated in any genome x refenrece comparisons from joint variant call datasets.

Because alignment depths don't scale with copy number (there is often NO correlation at all), it is impossible to identify illegals calls caused by repeated sequences using only the information contained in VCF files. SmartSNPs uses an external repeat finding tool and then applies three simple criteria to cull illegal SNP calls from VCF files:

1. The SNP should not be in a sequence that is repeated in the reference genome. SmartSNPs utilizes a novel de novo repeat-finding algorithm that outperforms related tools and creates lookup strings that record the copy number at each position in the reference genome. SmartSNPs cross-references SNP sites against the lookup strings and culls sites that occur in repeats (copy number > 1).

2. When variant calling in a haploid genome, the SNPs should not be heterozygous because this indicates that it resides in a repeat region in the query genome. Ideally variant callers should not call heterozygous SNPs when the "haploid" option is used but I guess some (e.g. the Genome Analysis Tool Kit - GATK) are not that smart! SmartSNPs culls sites where one or more samples possesses a heterozygous call. Some leeway is provided to allow for sequencing errors occasionally producing a small number of reads matching the reference allele.

3. Read depth should not be significantly lower or higher than the average depth for a single-copy sequence. Curiously, an occasional symptom of a SNP having been called in a repeated sequence is sequence coverage that is far lower than the average for single-scopy sequence (a function of the stochastic alignment of highly diverged sequences).

Running SmartSNPs

SmartSNPs is designed to run on individual SNP files owing to the difficulty of checking filtering 





Because alignment depths don't scale with copy number (there is often NO correlation at all), it is impossible to identify these illegals calls using only the information contained in VCF files. SmartSNPs uses external resources to identify filters SNP calls using the following criteria:

