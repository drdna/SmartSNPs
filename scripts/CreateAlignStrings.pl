#!/usr/bin/perl

## CreateAlignStrings.pl

## written by Mark L. Farman

## blasts a genome against itself and determines the copy number of each position in the genome


die "Usage: perl CreateAlignStrings.pl <ref_genome>\n" if @ARGV < 1;

use FetchGenomeOrig;


# Check fasta filename format

die "Reference genome filename can only contain one period before the final suffix\n" if $ARGV[0] =~ /\..+\./;


# create output directory name

($outdir = $ARGV[0]) =~ s/\.[A-Za-z]+?$//;


# create output directory

unless (-d "$outdir") {

  mkdir "$outdir" || die "can't make dir: $outdir, it already exists\n";

}

$GenomeLengthHashRef = REF_GENOME_LENGTHS($ARGV[0]);

$alignStringHashRef = ALIGN_STRINGS($GenomeLengthHashRef);

& SELFBLAST;

$alignStringHashRef = READ_BLAST($alignStringHashRef);

PRINT_STRINGS($alignStringHashRef);



## subroutines ##


# get lengths of genome chromosomes/contigs

sub REF_GENOME_LENGTHS {

  my $genome = $_[0];

  my $GenomeLengthHashRef = FetchGenomeOrig::getLens($genome);

}


# create the initial alignment strings showing 0 blast alignments at each site

sub ALIGN_STRINGS {

  my $GenomeLengthHashRef = $_[0];

  my %GenomeLengthHash = %$GenomeLengthHashRef;

  my %alignStringHash = undef;

  foreach my $Seq (keys %GenomeLengthHash) {

    $alignStringHash{$Seq} = "0" x $GenomeLengthHash{$Seq};

  }

  return \%alignStringHash;

}



# Blast genome against itself

sub SELFBLAST {

  $blastout = $outdir.".BLAST";

  $blastparams = "-query $ARGV[0] -subject $ARGV[0] -evalue 1e-20 -max_target_seqs 20000 -out $blastout -outfmt '6 qseqid sseqid qstart qend sstart send btop'";

  system("blastn $blastparams");

}



# read blast results
  
sub READ_BLAST {

  my $alignStringHashRef = $_[0];

  my %alignStringHash = %{$alignStringHashRef};

  my $AlignStringsOut = $outdir."_alignments";

  open(OUT, '>', "$outdir/$AlignStringsOut") || die "Can't open outfile\n";

  open(BLAST, "$blastout") || die "Can't open blast file\n";

  while(my $blast = <BLAST>) {

    ($qid, $sid, $qs, $qe, $ss, $se, $btop) = split(/\t/, $blast);

    substr($alignStringHash{$qid}, $qs-1, ($qe-$qs)+1) =~ tr/01/12/;

  }

  close BLAST;

  return \%alignStringHash

}


# print the alignment strings

sub PRINT_STRINGS {

  my $alignStringHashRef = $_[0];

  my %alignStringHash = %{$alignStringHashRef};

  foreach $record (sort {$a cmp $b} keys %alignStringHash) {

   next if $record eq '';

    print OUT "$record\t$alignStringHash{$record}\n"

  }

  close OUT;
   
}  
