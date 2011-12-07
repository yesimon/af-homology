#!/bin/perl

my $dna;

while (<STDIN>) {
	chomp;
	if ($_ =~ m/>/) {
		if (defined($dna)) {
			print "$dna\n";
		}
		$dna = "";
	} else {
		$dna = $dna . $_;
	}
}

print "$dna\n";
