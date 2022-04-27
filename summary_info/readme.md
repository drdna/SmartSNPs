**Summarizing SNP calls to identify problems that persist after filtering**

SNPs_summary.pl was used to read the filtered VCF files and collate the information in a format that is easily manually interrogated. The sample file provided here (SNPsummary.txt) reveals a number of problematic SNP calls that persist after filtering. A good example, is a massive run of SNPs on Chr3 in sample ERR2061623.

**Note:** It would be straightforward to automate this second phase of filtering but it is useful to perform this step manually at first to get a sense of the errors that still persist after filtering (note that these illegal calls are unavoidable and typically result from presence/absence polymorphisms that prevent even SmartSNPs from identifying repeat sequences):

	Usage: perl SNPs_summary.pl <VCF_directory>

A single output file is created (.summary.txt) which can be interrogated to make sure that SNPs do not occur in linkage blocks - these will typically show up as multiple consecutive SNPs found in the same group of strains (or a single strain). 
