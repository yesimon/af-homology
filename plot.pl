#!/usr/bin/env perl

my $target = $ARGV[0];

open ANS, "ans.txt" or die $!;

my $start;
my $end;
my $zeb_chr;
my $zeb_start;

while (<ANS>) {
	my @slice = split(/\s+/);
	my $cne = $slice[0];
	if (not ($cne eq $target)) {
		next;
	}

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

	$start = $ans_start - $zeb_start;
	$end = $start + $ans_end - $ans_start;
	last;
}

print "graphing $target...\n";

close ANS;

open (GNUPLOT, "|gnuplot");
print GNUPLOT "s(x) = x > $start && x < $end ? 8 : 0\n";
print GNUPLOT <<EOPLOT;
set terminal postscript enhanced "Courier" 14 linewidth 1 rounded

set style line 1 lc rgb '#000000' lt 1 lw 1
set style line 2 lc rgb '#ad0000' lt 1 lw 2
set output "$target.eps"
set size 1 ,1
set nokey
set xlabel "alignment"
set title "$target, $zeb_chr:$zeb_start (danRer5)"
set grid xtics ytics
plot "tmp.dat" w lines ls 1, s(x) ls 2
EOPLOT
close(GNUPLOT);
