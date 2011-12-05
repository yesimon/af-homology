#!/usr/bin/env perl
# hexymf_score.pl - Determines the score of a sequence
# usage: <score_result_file_name> <test_sequence>

my %score;

sub max {
    splice(@_, ($_[0] > $_[1]) ? 1 : 0, 1);
    return ($#_ == 0) ? $_[0] : max(@_);
}

sub count_mers {
	my ($str, $n) = @_;
	my %count = ();

	my $i;
	chomp($str);
	$str = uc($str);

	for ($i = 0; $i < length($str) - $n + 1; $i++) {
		my $key = substr($str, $i, $n);
		if (exists $count{$key}) {
			$count{$key} += 1;
		} else {
			$count{$key} = 1;
		}
	}

	return %count;
}

sub rc {
	my ($str) = @_;

	my $ans = "";
	chomp($str);
	my $rstr = reverse($str);
	@slice = split(//, $rstr);
	foreach $v (@slice) {
		my $r;

		if (uc($v) eq 'A') {
			$r = 'T';
		} elsif (uc($v) eq 'C') {
			$r = 'G';
		} elsif (uc($v) eq 'G') {
			$r = 'C';
		} elsif (uc($v) eq 'T') {
			$r = 'A';
		} else {
			#die "couldn't match nucleotide\n";
			$r = 'N';
		}

		$ans = $ans . $r;
	}

	return $ans;
}

sub ymf_load_training {
	my ($fname) = @_;

	%score = ();
	open SC, $fname or die $!;

	while (<SC>) {
		my @slice = split(/\s+/);
		my $key = $slice[0];

		if ($key eq "The") {
			next;
		} elsif ($key eq "0") {
			last;
		}

		my $zscore = $slice[2];
		$score{$key} = $zscore;
	}

	close SC;
}

sub ymf_train {
	my ($dna) = @_;
	my $dna_len = length($dna);
	system("echo \">\" > tmp.fa");
	open TMP, ">> tmp.fa" or die $!;
	print TMP "$dna\n";
	close TMP;
	system("./ymf_stats ymf.config $dna_len 6 human tmp.fa > /dev/null");
	system("rm tmp.fa");

	ymf_load_training("results");
	system("rm results");
}

sub ymf_score {
	my ($dna) = @_;
	my %count = count_mers($dna, 6);
	my $score = 0;

	foreach $val (sort {$count{$b} <=> $count{$a} } keys %count)
	{
		if (exists $score{$val}) {
			#print "$val\t$count{$val}\t$score{$val}\n";
			$score += $count{$val} * $score{$val};
		}
	}

	return $score;
}

while (<STDIN>) {
	my @slice = split(/\s+/);
	$cns = $slice[0];
	$chr = $slice[1];
	$start = $slice[2];
	$end = $slice[3];
	$zebcoord = $slice[4];
	$human_dna = $slice[5];
	$zebra_dna = $slice[7];

	@slice = split(/\=|\:|\-/, $zebcoord);
	my $zeb_chr = @slice[0];
	my $zeb_start = @slice[1];

	#print "human len is " . length($human_dna) . ", zebrafish len is " . length($zebra_dna) . "\n";

	ymf_train($human_dna);

	my $i;
	my $tlen = length($human_dna);
	for ($i = 0; $i < length($zebra_dna) - $tlen + 1; $i++) {
		my $dna = substr($zebra_dna, $i, $tlen);
#		my $score = max(ymf_score($dna), ymf_score(rc($dna)));
		my $score = ymf_score($dna . "N" . rc($dna));
		my $pstart = $zeb_start + $i;
		my $pend = $pstart + $tlen;
		print "$cns\t$chr\t$start\t$end\t$zeb_chr:$pstart-$pend\t$score\n";
	}
}

#foreach $value (sort {$count{$b} <=> $count{$a} } keys %count)
#{
#	print "$value\t$count{$value}\n";
#}
