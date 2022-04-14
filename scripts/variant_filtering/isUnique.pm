package isUnique;

my %alignHash = undef;

sub getAlignStrings {

  my ($alignFile) = @_;

  open(ALIGNMENT, "$alignFile") || die "Can't open align file $alignFile\n";

  while(my $L = <ALIGNMENT>) {

    chomp($L);

    next if $L eq '';

    my($contig, $alignString) = split(/\t/, $L);

    $contig =~ s/.+?(\d+)$/$1/;

    $alignHash{$contig} = $alignString;

  }

  close ALIGNMENT;

  return \%alignHash

}

sub unique {

  my ($alignFile, $contig, $position) = @_;

  open(ALIGNMENT, "$alignFile") || die "Can't open align file $alignFile\n";

  while(my $L = <ALIGNMENT>) {

    chomp($L);

    next if $L eq '';

    my($contig, $alignString) = split(/\t/, $L);

    $alignHash{$contig} = $alignString;

  }

  if(substr($alignHash{$contig}, $position, 1) == 1) {

    return 'y'

  }

  else {

    return 'n'

  }

  close ALIGNMENT

}

1; 
