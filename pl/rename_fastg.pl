#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);


if ($#ARGV == -1) {
        print "usage: $0  <input fastg> <input fasta> <output fastg>\n";
        exit 1;
}

my $ifn_fastg = $ARGV[0];
my $ifn_fasta = $ARGV[1];
my $ofn = $ARGV[2];

####################################################################################
# read fasta
####################################################################################

my %fasta;
my %count2contig;

print "reading fasta file: $ifn_fasta\n";
open(IN, $ifn_fasta) || die "cannot open $ifn_fasta";

my $count = 0;
my $seq = "";
my $contig = "";
while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) ne ">") {
	$seq .= $line;
    } else {
	if ($contig ne "") {
	    $fasta{$contig} = {};
	    $fasta{$contig}->{seq} = $seq;
	    $count2contig{$count} = $contig;
	    $count++;
	}
	my @f = split(" ", substr($line,1));
	$contig = $f[0];
	$seq = "";
    }
}
if ($contig ne "") {
    $fasta{$contig} = {};
    $fasta{$contig}->{seq} = $seq;
    $count2contig{$count} = $contig;
}
close(IN);

print "number of contigs: ", scalar keys %fasta, "\n";

####################################################################################
# traverse fastg, first time
####################################################################################

my %map;

print "traversing fastg file to build map: $ifn_fastg\n";
open(IN, $ifn_fastg) || die "cannot open $ifn_fastg";

$count = 0;
$seq = "";
$contig = "";
while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) ne ">") {
	$seq .= $line;
    } else {
	if ($contig ne "" && !($contig =~ /'$/)) {
	    defined($count2contig{$count}) or die $count;
	    my $fcontig = $count2contig{$count};
	    $fasta{$fcontig}->{seq} eq $seq or
		die "sequence mismatch, make sure fasta and fastg are compatible";
	    $map{$contig} = $fcontig;
	    $count++;
	}
	chop($line);
	my @f = split(":", substr($line,1));
	$contig = $f[0];
	$seq = "";
    }
}
if ($contig ne "" && !($contig =~ /'$/)) {
    defined($count2contig{$count}) or die $count;
    my $fcontig = $count2contig{$count};
    $fasta{$fcontig}->{seq} eq $seq or
	die "sequence mismatch, make sure fasta and fastg are compatible";
    $map{$contig} = $fcontig;
    $count++;
}
close(IN);

####################################################################################
# traverse fastg, second time
####################################################################################

print "re-traversing fastg file: $ifn_fastg\n";
open(IN, $ifn_fastg) || die "cannot open $ifn_fastg";

print "writing fastg file: $ofn\n";
open(OUT, ">", $ofn) || die "cannot create $ofn";

while (my $line = <IN>) {
    chomp($line);
    if (substr($line, 0, 1) ne ">") {
	print OUT $line, "\n";
    } else {
	chop($line);
	my @f = split(":", substr($line,1));

	# sanity checks
	scalar @f > 0 or die "cannot find ':' in line: $line"; 
	scalar @f <=2 or die "expecting single ':' in line: $line";

	# handle source contig
	my $src_contig = $f[0];
	my $hline = ">".rename_contig($src_contig);

	if (scalar @f > 1) {
	    $hline .= ":";
	    my @tgt_contigs = split(",", $f[1]);
	    my $N = scalar @tgt_contigs;
	    for (my $i=0; $i<$N; $i++) {
		$hline .= "," if ($i > 0);
		my $cc = rename_contig($tgt_contigs[$i]);
		$hline .= $cc;
	    }
	}
	print OUT $hline, ";\n";
    }
}
close(IN);
close(OUT);

######################################################################################################
# Subroutines
######################################################################################################

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

sub rename_contig
{
    my ($contig) = @_;
    my $is_reverse = $contig =~ /'$/;
    $contig = substr($contig, 0, length($contig)-1) if ($is_reverse);
    defined($map{$contig}) or die "contig $contig not defined in map";
    
    my $result = $map{$contig};
    $result .= "'" if ($is_reverse);
    print $contig, " -> ", $result, "\n";
    return $result;
}
