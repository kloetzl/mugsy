#!/usr/bin/perl
#Returns list of blocks that contain all sequences in the set seqid1...seqidn
#./mafgrep.pl seqid1 seqid2 ... seqidn < out.maf

use strict;

my $format='maf';#or tab

my %grepids = map { $_, 1 } @ARGV;
print STDERR "Looking for ",scalar(keys %grepids),"\n";
my $currscore;
my $currorient;
my $blockorient;
my @allblocks;
my $block = [];
while(my $line=<STDIN>){
    if($line =~ /^a\s+score=(\S+)/){
	$currscore=$1;
	push @allblocks,$block;
	$block=[];
    }
    elsif($line =~ /^s/){
	my @elts = split(/\s+/,$line);
	#0-score,1-blockorient,2-accession,3-start,4-end
	push @$block,[$currscore,$currorient,$elts[1],$elts[2],$elts[2]+$elts[3],$elts[3],$elts[4],$line];
    }
}
print STDERR "Parsed ",scalar(@allblocks)," blocks\n";
print "##maf version=12\n";
push @allblocks,$block;
foreach my $blocks (@allblocks){
    #Lookup of all seqs in the block
    my %seqs = map {$_->[2], 1} @$blocks;
    #
    my %results = map { $_, $grepids{$_} } grep { not exists $seqs{$_} } keys %grepids;
    #print STDERR "Seqs ",join(',',sort keys %seqs)," ",scalar(@$block),"\n";
    #print STDERR "Results ",join(',',sort keys %results),"\n";
    #print STDERR "Grep ",join(',',sort keys %grepids),"\n";	#join(' ',keys %seqs)," | ",join(' ',keys %grepids),"\n";
    if(scalar(keys %results)==0){
	if($format eq 'maf'){
	    print "a score=$blocks->[0]->[0]\n";
	}
	foreach my $bl (@$blocks){
	    if($format eq 'maf'){
		if(exists $grepids{$bl->[2]}){
		    print "$bl->[7]";
		}
	    }
	    else{
		print "$bl->[2]\t$bl->[3]\t$bl->[4]\t$bl->[5]\t$bl->[6]\n";
	    }
	}
	print "\n";
    }
}
