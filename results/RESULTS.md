**SNP Filtering Pipeline Overview:**

1.	Use **SmartSNPsV2.pl** to filter out illegal SNPs

Briefly, SmartSnpsV2.pl operates on INDIVIDUAL vcf files (not joint calls) as follows:

a)	Reads an alignment file that lists whether a nucleotide position is unique (value = 1) or repeated (value = 2) and hashes each alignment
b)	Compares each SNP position with the corresponding alignment string and rejects if it occurs in a chromosome position where the copy number is > 1 (i.e. a repeated region)
c)	Examines the genotype information for each SNP and rejects is the Alt:Ref ratio is not greater a user-specified value (I suggest using 10). This builds in a “cushion” to allow for very occasional sequencing errors.
d)	Determine high quality read coverage across the SNP position (alt + ref) and rejects the SNP if the number of reads is less than a user-defined cut-off (I suggest 20 but this will depend on overall read coverage in a given dataset. Good SNPs that have been rejected can be added back in using a SNP_panning.pl script. However, experience tells us that most rejected sites are problematic anyway).

  perl SmartSNPs/scripts/variant_filtering/SmartSNPsV2.pl <VCF> \ BDIR_align/B71v2.B71v2_alignments 20 20 1


Two output files are created: 

a)	a VCF file (.SSfilter.vcf) with PASS/FAIL filter results
b)	a .log file that list the criterion used to reject each SNP.


2.	Use **SNPsummary.pl** to read the filtered VCF files and collate the information in a format that is easily manually interrogated. Will be straightforward to automate this second phase of filtering but it is useful to perform this step manually at first to get a sense of the errors that still remain after filtering (note that these illegal calls are unavoidable and typically result from presence/absence polymorphisms that prevent even SmartSNPs from identifying repeat sequences):

  perl SNPsummary.pl <VCF_directory>

A single output file is created (summary.txt) which can be interrogated to make sure that SNPs do not occur in linkage blocks (these will typically show up as multiple consecutive SNPs found in the same group of strains (or a single strain). 

3.	Iterate through filtered VCF files using different values for Alt:Ref ratios/low coverage cut-off (10 to 50, step 5), capture summary metrics and write to output files:

lowCoverage cutoff:

  for f in `ls BGL_ZMB_VCFs/*snps.vcf`; do for g in 10 15 20 25 30 35 40 45 50; do suff=${f/*\//}; pref=${suff/\.*/}; perl SmartSNPs/scripts/variant_filtering/SmartSNPsV2.pl $f BDIR_align/B71v2.B71v2_alignments 30 $g 1 | grep Allowed | awk -v var="$g" -v pref="$pref" -F '; |: ' '$3 ~ /Allowed/ {OFS="\t"} {print var, $4, $6, $8, $10}' >> low.txt; done; done

Alt:Ref ratio:

  for f in `ls BGL_ZMB_VCFs/*snps.vcf`; do for g in 10 15 20 25 30 35 40 45 50; do suff=${f/*\//}; pref=${suff/_*/}; perl SmartSNPs/scripts/variant_filtering/SmartSNPsV2.pl $f BDIR_align/B71v2.B71v2_alignments $g 20 1 | grep Allowed | awk -v var="$g" -v pref="$pref" -F '; |: ' '$3 ~ /Allowed/ {OFS="\t"} {print pref, var, $4, $6, $8, $10}' >> heterozygous.txt; done; done

4.	Use **SNPfilter.R** to plot summary statistics and look for shoulders as indicators of best values to select for final filtering.


Notes:

While, in principle, it would be possible to modify the SmartSNPsV2.pl script to run on merged VCF files, it is virtually impossible to run manual checks for anomalies on the filtered data.
