# SNP filtering pipeline

# Mark L Farman, 04/14/22


refgenome=$1

VCFdir=$2

ratio=$3

depth=$4

cutoff=$5

perl SmartSNPs/scripts/CreateAlignStrings.pl $refgenome;

alignfile=${refgenome/\.*/}/${refgenome/\.*/}_alignments

for f in `ls $VCFdir/*vcf`; do echo $f; perl SmartSNPs/scripts/SmartSNPs.pl $f $alignfile $ratio $depth $cutoff; done




