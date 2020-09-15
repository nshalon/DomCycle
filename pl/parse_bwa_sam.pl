#!/usr/bin/env perl

use strict;
use POSIX;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <ifn> <ofn> <ofn stats>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];
my $ofn_stats = $ARGV[2];

######################################################################################################
# output syntax

# id: sequence id
# contig, start_coord, end_coord, strand: reference bounding region
# edit_dist: edit distance
# score: map quality
# match_length: total matching nt length
# cigar: cigar string
# subsitute: substitute nt for nt
#            ;(reference_coord,sequence_coord,reference_nt,sequence_nt;)*
# insert: insert nts before ref coord
#         ;(reference_coord,seqeunce_coord,nts;)*
# delete: delete nts on coord
#         ;(reference_coord,length;)*
# clip: clip start and end of read sequence (soft clip)
#       start;end of clip in sequence coords
# sequence: read sequence
######################################################################################################

print STDERR "traversing read file: $ifn\n";
open(IN, $ifn) || die $ifn;

print STDERR "generating file: $ofn\n";
open(OUT, ">", $ofn) || die $ofn;
print OUT "id\tcontig\tcoord\tback_coord\tstrand\tedit_dist\tscore\tmatch_length\tcigar\tsubstitute\tinsert\tdelete\tclip\tunique\tsequence\n";

my $count = 0;

my %stats;
$stats{ok} = 0;
$stats{multi_segment} = 0;
$stats{unmapped} = 0;
$stats{duplicate} = 0;

while (my $line = <IN>) {
    chomp($line);

    # skip comments
    next if (substr($line, 0, 1) eq "@");

    $count++;

    print STDERR "line: $count\n" if ($count % 100000 == 0);
    my @f = split(/\s+/, $line);
    my $id = $f[0];
    my $flag = $f[1];
    my $contig = $f[2];
    my $coord = $f[3];
    my $score = $f[4];
    my $cigar = $f[5];
    my $seq = $f[9];

    my $strand = ($flag & 16) ? -1 : 1;
    my $supp_alignment = ($flag & 2048);

    # skip all kinds of bad reads
    my $multi_segments = ($flag & 1);
    my $unmapped = ($flag & 4);
    my $duplicate = ($flag & 1024);
    if ($multi_segments || $unmapped || $duplicate) {
	if ($multi_segments) {
	    $stats{multi_segment}++;
	} elsif ($unmapped) {
	    $stats{unmapped}++;
	} else {
	    $stats{duplicate}++;
	}
	next;
    }
    $stats{ok}++;

    print STDERR "Warning, unexpected flag: $flag\n" if ($flag & ~(16 | 2048));

    my @read2ref = (0) x length($seq);
    my @match2read = (0) x length($seq);

    # get edit distance
    my $dist = "";
    my $md = "";
    my $unique = 1;
    for (my $i = 11; $i < @f; $i++) {
	$dist = $f[$i] if substr($f[$i], 0, 2) eq "NM";
	$md = $f[$i] if substr($f[$i], 0, 2) eq "MD";
	$unique = 0 if substr($f[$i], 0, 2) eq "XA";
    }
    $dist =~ s/NM:i://;
    $dist ne "" or die "NM field not found";

    # insertion and deletion
    my %hdelete;
    my $delete_str = ";";
    my $insert_str = ";";

    # first go over cigar and setup maps between spaces
    my $rcoord = $coord;
    my $scoord = 0;
    my $mcoord = 0;
    my @fcigar = $cigar =~ m/\d+\w/g;
    for (my $i=0; $i<@fcigar; $i++) {
	my ($len, $type) = $fcigar[$i] =~ m/(\d+)(\w)/g;
        # print "CIGAR m:$mcoord, s:$scoord, r:$rcoord\n";

	if ($type eq "D") {
	    $hdelete{$rcoord} = $len;
	    # print "D:",$rcoord,",",$len,"\n";
	    $delete_str .= $rcoord.",".$len.";";
	    $rcoord += $len;
	}

	if ($type eq "I") {
	    my $seq_nts = substr($seq, $scoord, $len);
	    $insert_str .= $rcoord.",".($scoord+1).",".$seq_nts.";";
	}

	$scoord += $len if ($type eq "S" || $type eq "I");

	if ($type eq "M") {
	    for (my $j = 0; $j < $len; $j++) {
		$match2read[$mcoord+$j] = $scoord + $j;
		$read2ref[$scoord+$j] = $rcoord + $j;
	    }
	    $mcoord += $len;
	    $scoord += $len;
	    $rcoord += $len;
	}
    }
    # total match length
    my $mlength = $mcoord;

    # back_coord/front_coord are the first/last coord of the map pointing towards the sequenced molecule
    my $front_coord = ($strand == 1 ? $rcoord-1 : $coord);
    my $back_coord = ($strand == 1 ? $coord : $rcoord-1);

    # substitutions are specifed in the MD field, in match space
    my $sub_str = ";";
    $mcoord = 0;
    $rcoord = $coord;
    ($md ne "" and substr($md, 0, 5) eq "MD:Z:") or die "MD not found";
    $md =~ s/MD:Z://;
    while ($md =~ /([0-9]+|[A-Z]|\^[A-Z]+)/g) {
	my $mfield = $1;
        # print "FIELD: ", $mfield, "\n";

        # print "MD m:$mcoord, s:$scoord, r:$rcoord\n";

	if ($mfield =~ /[0-9]+/) {
            # match
	    $mcoord += $mfield;
	    $rcoord += $mfield;
	} elsif ($mfield =~ /\^([A-Z]+)/) {
	    # deletion
	    my $length = length($mfield)-1;
	    # print "Dx:",$rcoord,",",$length,"\n";

	    (defined($hdelete{$rcoord}) and $hdelete{$rcoord} == $length) or die "cigar and MD field do not match on delete";
	    $rcoord += $length;
	} elsif ($mfield =~ /([A-Z])/) {
	    # mismatch
	    my $ref_nt = $mfield;
	    my $scoord = $match2read[$mcoord];
	    my $seq_nt = substr($seq, $scoord, 1);

	    $sub_str .= $rcoord.",".($scoord+1).",".$ref_nt.",".$seq_nt.";";
	    $mcoord++;
	    $rcoord++;
	} else {
	    die "failed parsing MD field";
	}
    }

    # soft clips
    my ($len0, $type0) = $fcigar[0] =~ m/(\d+)(\w)/g;
    my $clip_start = ($type0 eq "S") ? ($len0+1) : 1;

    my ($lenN, $typeN) = $fcigar[$#fcigar] =~ m/(\d+)(\w)/g;
    my $clip_end = ($typeN eq "S") ? length($seq)-$lenN : length($seq);

    my $clip_str = $clip_start.";".$clip_end;

    my $unique_str = $unique ? "T" : "F";
    print OUT "$id\t$contig\t$front_coord\t$back_coord\t$strand\t$dist\t$score\t$mlength\t$cigar\t$sub_str\t$insert_str\t$delete_str\t$clip_str\t$unique_str\t$seq\n";

}
close(IN);
close(OUT);

print_hash($ofn_stats, %stats);

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

