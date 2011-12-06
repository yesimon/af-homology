#!/usr/bin/env perl
# hexdiff_score.pl - Determines the score of a sequence
# usage: <score_result_file_name> <test_sequence>

#use List::Util 'shuffle';

my %score;

sub max {
    splice(@_, ($_[0] > $_[1]) ? 1 : 0, 1);
    return ($#_ == 0) ? $_[0] : max(@_);
}

sub shuffleStr {
    my $len = length $_[0];
    my ($tmp, $n);
    $n = $_+rand($len-$_)
    , $tmp = substr( $_[0], $_, 1)
    , substr( $_[0], $_, 1) = substr( $_[0], $n , 1)
    , substr( $_[0], $n , 1) = $tmp
        for 0 .. $len;
    $_[0];
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

sub hexdiff_train {
	my ($dna) = @_;

	%score = ();

	my %test = count_mers($dna, 6);
	my %background = count_mers(shuffleStr($dna), 6);

	foreach $val (sort {$test{$b} <=> $test{$a} } keys %test) {
		#print "$val\t$test{$val}\t$background{$val}\n";
		my $q = exists($background{$val}) ? $test{$val} / $background{$val} : $test{$val};

		if ($q > 1) {
			#print "$val\t$q\n";
			$score{$val} = $q;
		}
	}
}

sub ymf_score {
	my ($dna) = @_;
	my %count = count_mers($dna, 6);
	my $score = 0;

#	foreach $val (sort {$count{$b} <=> $count{$a} } keys %count)
	foreach $val (keys %count)
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
	my $cns = $slice[0];
	my $chr = $slice[1];
	my $start = $slice[2];
	my $end = $slice[3];
	my $zebcoord = $slice[4];
	my $human_dna = $slice[5];
	my $zebra_dna = $slice[7];

	@slice = split(/\=|\:|\-/, $zebcoord);
	my $zeb_chr = @slice[0];
	my $zeb_start = @slice[1];

	hexdiff_train($human_dna);

	my $i;
	my $tlen = length($human_dna);
	for ($i = 0; $i < length($zebra_dna) - $tlen + 1; $i++) {
		my $dna = substr($zebra_dna, $i, $tlen);
		my $score = max(ymf_score($dna), ymf_score(rc($dna)));
		#my $score = ymf_score($dna . "N" . rc($dna));
		my $pstart = $zeb_start + $i;
		my $pend = $pstart + $tlen;
		print "$cns\t$chr\t$start\t$end\t$zeb_chr:$pstart-$pend\t$score\n";
	}
}

