#!/usr/bin/perl -w

use strict;
use warnings;

my %aminoacids = (
	'ALA', 'A',
	'ARG', 'R',
	'ASN', 'N',
	'ASP', 'D',
	'CYS', 'C',
	'GLN', 'Q',
	'GLU', 'E',
	'GLY', 'G',
	'HIS', 'H',
	'ILE', 'I',
	'LEU', 'L',
	'LYS', 'K',
	'MET', 'M',
	'PHE', 'F',
	'PRO', 'P',
	'SER', 'S',
	'THR', 'T',
	'TRP', 'W',
	'TYR', 'Y',
	'VAL', 'V' );
my $seq;
my $temp;
my $infile = $ARGV[0] || die "# usage: $0 <promoters file>\n";
open(SEQ, $infile) || die "# cannot open input $infile : $!\n";
while(<SEQ>){
	if (/^(\w+)/) {
		$temp = $aminoacids{$1};
		$seq = $seq.$temp;
	}
}
close(SEQ);
print "$seq\n";
