#!/usr/bin/perl

use strict;
use warnings;

# Argument checking
if (scalar @ARGV < 2) {
  print STDERR "Usage:\n removeBad-dev.pl inputDEV chr position\n";
  print STDERR "  Removes part of dev that comes from a single call with chr/position into\na new devfiles with a .trimmed.dev suffix\n";
  print STDERR "\n";
  exit(1);
}

my $devFile = $ARGV[0];
my $chrOfCallToRemove = $ARGV[1];
my $posOfCallToRemove = $ARGV[2];

# HACK, assumes chromosome names start with chr
# Adjust accordingly
my $startVCFLine = "chr";

my $printRecord = 1;  # continue to print record

#Read in dev file, and print out most items
open(DEV, $devFile) or die $!;
open(NEWDEV, ">$devFile.trimmed.dev") or die $!;

while(<DEV>) {
  my $line = $_;
  if($line =~ /^$startVCFLine/) {
      if($printRecord==1) {
	  my @columns = split(/\t/,$line);
	  my ($chr,$pos) = ($columns[0],$columns[1]);
	  if($chr eq $chrOfCallToRemove && $pos==$posOfCallToRemove) {
	      $printRecord = 0;
	  }
      } else {
	  $printRecord = 1;
      }
  }
  print NEWDEV $line if $printRecord;
}
