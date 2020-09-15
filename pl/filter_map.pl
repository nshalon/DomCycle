#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <min quality score> <min length> <min edit distance> <ofn> <ofn stats>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $min_score = $ARGV[1];
my $min_length = $ARGV[2];
my $min_distance = $ARGV[3];
my $ofn = $ARGV[4];
my $ofn_stats = $ARGV[5];

print "min score: $min_score\n";
print "min length: $min_length\n";
print "min distance: $min_distance\n";

# stats
my $count = 0;
my $kept = 0;

print STDERR "traversing read file: $ifn\n";
open(IN, $ifn) || die $ifn;
my $header = <IN>;
chomp($header);
my %h = parse_header($header);

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;

my @fh = split("\t", $header);
$header = join("\t", @fh[0...$#fh-1]);
print OUT $header, "\n";

my %stats;
$stats{ok} = 0;
$stats{low_score} = 0;
$stats{high_edit_dist} = 0;
$stats{short_match_length} = 0;

while (my $line = <IN>) {
    chomp($line);

    $count++;
    print $count, "\n" if ($count % 1000000 == 0);

    my @f = split("\t", $line);
    $line = join("\t", @f[0...$#f-1]);

    # low quality (either non-unique or many too many mismatches in alignment)
    if ($f[$h{score}] < $min_score && $f[$h{unique}] eq "T") {
	$stats{low_score}++;
	next;
    }

    # high edit distance too high
    if ($f[$h{edit_dist}] > $min_distance) {
	$stats{high_edit_dist}++;
	next;
    }

    # matched sequence is too short
    if ($f[$h{match_length}] < $min_length) {
	$stats{short_match_length}++;
	next;
    }
    $stats{ok}++;
    print OUT $line, "\n";
}
close(IN);
close(OUT);

print_hash($ofn_stats, %stats);

######################################################################################################
# Subroutines
######################################################################################################

sub print_hash
{
    my ($ofn, %h) = @_;

    print STDERR "generating file: $ofn\n";
    open (OUT, ">", $ofn) || die $ofn;

    my $first = 1;
    foreach my $key (keys %h) {
	if ($first) {
	    print OUT $key;
	    $first = 0;
	} else {
	    print OUT "\t", $key;
	}
    }
    print OUT "\n";
    $first = 1;
    foreach my $key (keys %h) {
	if ($first) {
	    print OUT $h{$key};
	    $first = 0;
	} else {
	    print OUT "\t", $h{$key};
	}
    }
    print OUT "\n";
    close(OUT);
}

sub parse_header
{
	my ($header) = @_;
	chomp($header);
	my @f = split("\t", $header);
	my %result;
	for (my $i = 0; $i <= $#f; $i++) {
		$result{$f[$i]} = $i;
	}
	return %result;
}
