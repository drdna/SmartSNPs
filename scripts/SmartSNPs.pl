#!/usr/bin/perl

## SmartSNPs.pl: a program to filter illegal SNPs calls in VCF files

## written by Mark L. Farman


## Usage warnings

die "Usage: perl SmartSNPs.pl <VCF> <alignfile> <AltRefRatio> <ReadDepth> [<cutoff>]\n" if @ARGV < 4;


## load module

use isUnique;


## Read arguments

($VCF, $alignFile, $ratio, $depth, $cutoff) = @ARGV if @ARGV == 5;

($VCF, $alignFile, $ratio, $depth) = @ARGV if @ARGV == 4;

$cutoff = 1 if $cutoff eq '';



## read in alignment strings that record copy number at each site

& READ_ALIGNSTRINGS;

## open VCF file and read sample info line

open(VCF, $VCF);

print "$VCF\n";		# prints a VCF file ID header for each set of retained sites

$data = 'no';

while($V = <VCF>) {

  chomp($V);

  if($V =~ /^\#CHROM/) {

    $data = 'yes';	# flag indicating data lines are starting

    next

  }

#  print "$V\n" if $data eq 'no';	# uncomment this line if the complete VCF header is required in the output

  if($data eq 'yes') {

    $numRecords ++;			# count total variants in file

    ($ID, $pos, $varid, $ref, $alt, $qual, $filter, $info, $format, $genotypeInfo) = split(/\t/, $V, 10);

    $ID =~ s/.+?(\d+)/$1/;
   
    & REPEATS;

    & HETEROZYGOTES if $cutoff < 2;

    & LOW_COVERAGE;

#    print "$V\n" if $allowed eq 'yes';
  
  }

}



## SUBROUTINES ##


## SUBROUTINES

sub READ_ALIGNSTRINGS {

  # use isUnique module to read in the alignments strings

  $AlignHashRef = isUnique::getAlignStrings($alignFile);

  %AlignHash = %$AlignHashRef;

}


sub REPEATS {

  # look at SNP site in corresponding alignment string
  # reject site if # of alignments is greater than specified cutoff
  # increment the repeated sites counter

  $allowed = 'yes';

  if(substr($AlignHash{$ID}, $pos-1, 1) > $cutoff) {

    $allowed = 'no';

    $substr = substr($AlignHash{$ID}, $pos-1, 1);

    $repeat ++;

    next

  }

}


sub HETEROZYGOTES {

  # check for heterozygous genotype calls
  # reject site if number of alt allele calls is less than 50 times ref allele calls 
  # cutoff needs to be adjusted based on overall read depth

  $allowed = 'yes';

  @genotypeInfo = split(/:/, $genotypeInfo);

  ($ref, $alt) = split(/,/, $genotypeInfo[1]);

  if($genotypeInfo[0] == 1 && $ref > 0) {

    unless($alt/$ref >= $ratio) {

      print "$V\n";	# uncomment this line to view rejected VCF lines

      $heterozygote ++;

      $allowed = 'no';

      next

    }

  }

}


sub LOW_COVERAGE {

  # reject line if read depth used to infer genotype is < 30
  # adjust according to average read depth

  $allowed = 'yes';

  if($ref+$alt < $depth) {

    $lowCoverage ++;

    $allowed = 'no';

    next

  }

}

## calculate summary statistics

$allowed = $numRecords - $repeat - $heterozygote - $lowCoverage;


## print summary statistics

print "#NumRecords: $numRecords; Allowed: $allowed; Repeated: $repeat; Non-repeat heterozygotes: $heterozygote; Low coverage: $lowCoverage\n" if $cutoff == 1;

print "#NumRecords: $numRecords; Allowed: $allowed; Repeated: $repeat; Low coverage: $lowCoverage\n" if $cutoff > 1;

close VCF;
