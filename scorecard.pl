#!/usr/bin/env perl 

sub calc_overlap {
        my ($a_start, $a_end, $b_start, $b_end) = @_;

        if ($a_start >= $b_start && $a_start < $b_end) {
                my $tmp1 = $b_end - $a_start;
                my $tmp2 = $a_end - $a_start;

                return ($tmp1, $tmp2)[$tmp1 > $tmp2];
        } elsif ($a_start < $b_start && $a_end > $b_start) {
                my $tmp1 = $a_end - $b_start;
                my $tmp2 = $b_end - $b_start;

                return ($tmp1, $tmp2)[$tmp1 > $tmp2];
        }

        return 0;
}

open ANS, "ans.txt" or die $!;

my %astart;
my %aend;
my %alen;

while (<ANS>) {
	my @slice = split(/\s+/);
	my $cne = $slice[0];
	my $zebcoord = $slice[1];
	my $human_chr = $slice[2];
	my $human_start = $slice[3];
	my $human_end = $slice[4];
	my $anscoord = $slice[5];

	@slice = split(/\=|\:|\-/, $zebcoord);
	my $zeb_chr = @slice[0];
	my $zeb_start = @slice[1];
	my $zeb_end = @slice[2];

	@slice = split(/\=|\:|\-/, $anscoord);
	my $ans_chr = @slice[0];
	my $ans_start = @slice[1];
	my $ans_end = @slice[2];

	my $start = $ans_start - $zeb_start;
	my $end = $start + $ans_end - $ans_start;

	$astart{$cne} = $start;
	$aend{$cne} = $end;
	$alen{$cne} = $human_end - $human_start;
}

my $i = 0;
my $max_val = 0;
my $max_pos = 0;
my $last_cne;

my $wins = 0;
my $loses = 0;

my $win_names = "";

while (<STDIN>) {
	my @slice = split(/\s+/);
	my $cne = $slice[0];
	my $score = $slice[1];

	if (!defined($last_cne)) {
		$last_cne = $cne;
	}

	if (not ($cne eq $last_cne)) {
		if (calc_overlap($astart{$last_cne}, $aend{$last_cne},
				 $max_pos, $max_pos + $alen{$last_cne})) {
			$wins += 1;
			$win_names = $win_names . $last_cne . " ";
		} else {
			$loses += 1;
		}

		$i = 0;
		$max_val = 0;
		$max_pos = 0;
		$last_cne = $cne;
	}

	if ($score > $max_val) {
		$max_val = $score;
		$max_pos = $i;
	}

	$i += 1;
}

if (calc_overlap($astart{$last_cne}, $aend{$last_cne},
		 $max_pos, $max_pos + $alen{$last_cne})) {
	$wins += 1;
	$win_names = $win_names . $last_cne;
} else {
	$loses += 1;
}

my $total = $wins + $loses;

print "SCORE CARD: Got $wins/$total enhancers correct\n";

print "They are as follows: $win_names\n";
