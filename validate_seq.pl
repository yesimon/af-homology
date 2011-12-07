#!/usr/bin/perl

while (<STDIN>) {
	my @slice = split(/\s+/);
	my $cns = $slice[0];
	my $chr = $slice[1];
	my $start = $slice[2];
	my $end = $slice[3];
	my $zebcoord = $slice[4];
	my $human_dna = $slice[5];
	my $human_region = $slice[6];
	my $zebra_dna = $slice[7];

	@slice = split(/\=|\:|\-/, $zebcoord);
	my $zeb_chr = @slice[0];
	my $zeb_start = @slice[1];
	my $zeb_end = @slice[2];

	my $human_len = $end - $start;
	my $zeb_len = $zeb_end - $zeb_start;

	print "$human_len $zeb_len $zeb_len. Is " . length($human_dna) . " " . length($human_region) . " " . length($zebra_dna) . "\n";


	if ($human_len != length($human_dna) ||
	    $zeb_len != length($zebra_dna)) {
		die "Should be human enhanzer $human_len, human region $zeb_len, zebrafish region $zeb_len. Is " . length($human_dna) . " " . length($human_region) . " " . length($zebra_dna) . "\n";
	}
}

print "All tests passed.\n";
