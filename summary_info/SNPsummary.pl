#!/usr/bin/perl

die "Usage: SNPsummary.pl <VCF-dir>\n" if @ARGV < 1;

die "Program expects a DIRECTORY as argument\n" unless -d $ARGV[0];


## create output file

($outPrefix = $ARGV[0]) =~ s/.+\///;

$outPrefix =~ s/\/$//;

$outfile = $outPrefix.".summary.txt";

open(SUMMARY, '>', "$ARGV[0]/$outfile");

print SUMMARY join ("\t", ('contig', 'position', '#samples', 'samples_with_snp')), "\n";


## start reading VCF files

opendir(DIR, $ARGV[0]) || die "Can't open directory\n";

@filesList = readdir(DIR);

foreach $file (@filesList) {

  next unless $file =~ /SSfilter/;

  open(F, "$ARGV[0]/$file")|| die "Can't open file\n";

  while($L=<F>) {

    next if $L =~ /^#/;

    next if $L =~ /FAIL/;

    chomp($L);

    @L = split(/\t/, $L);

    if ($file =~ /(ERR206\d+|ERR2188\d+)/ && $ID !~ /(25|26)_.+vcf/) {		# skip over non-wheat blast datasets

      $ERR = $1;

      push @{$Hash{$L[0]}{$L[1]}}, $ERR;

    }

    elsif($file =~ /(ERR552\d+)/) {

      $ERR = $1;

      push @{$Hash{$L[0]}{$L[1]}}, $ERR;

    }

  }

}

  
foreach $chr (sort {$a cmp $b} keys %Hash) {
 
  foreach $pos (sort {$a <=> $b} keys %{$Hash{$chr}}) {

    $BGLsnps++ if $pop eq 'BGL';

    $ZMBsnps++ if $pop eq 'ZMB';

    $arrayLen = @{$Hash{$chr}{$pos}};

    print SUMMARY "$chr\t$pos\t$arrayLen\t@{$Hash{$chr}{$pos}}\n";

  }

}

print SUMMARY "BGL SNPs: $BGLsnps\n" if $pop eq 'BGL';         # uncomment for total SNPs summaries

print SUMMARY "ZMB SNPs: $ZMBsnps\n" if $pop eq 'ZMB';         # uncomment for total SNPs summaries
