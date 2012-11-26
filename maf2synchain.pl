#!/usr/bin/perl
#Convert MAF to FASTA
#Optionally only convert blocks that contain label 
#./maf2fasta.pl [label] < maf > fasta

use strict;

my $anchors = {};
my $seq2anchors = {};
my $seq2index = {};
my $genome2index = {};

my $anchornum=-1;
while(my $line=<STDIN>){
    if($line =~ /a\s+score=([\d\.\-]+)/){
	$anchornum++;
    }
    elsif($line =~ /s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\S+)/){
	my $accession = $1; #Must be formated as Genome.Sequence
	my $start = $2;
	my $len = $3;
	my $orientation = $4;
	my $end;
	if($orientation eq '-'){
	    $end = $start-$len-1;
	}
	else{
	    $end = $start+$len;
	}
	my $seqlength = $5;
	my $sequence;
	my $genome;
	if($accession =~ /([^\.]+)\.(\S+)/){
	    $genome=$1;
	    $sequence=$2;
	}
	else{
	    die "Accession not in Genome.Sequence format";
	}
	#Store index for this accession if first time we've seen it
	if(!exists $seq2index->{$accession}){
	    $seq2index->{$accession} = scalar(keys %$seq2index);
	}
	if(!exists $genome2index->{$genome}){
	    $genome2index->{$genome} = scalar(keys %$genome2index);
	}
	$anchors->{$anchornum}->{$accession}->{'gidx'} = $genome2index->{$genome};
	$anchors->{$anchornum}->{$accession}->{'sidx'} = $seq2index->{$accession};
	$anchors->{$anchornum}->{$accession}->{'start'} = ($start<$end ? $start:$end);
	$anchors->{$anchornum}->{$accession}->{'end'} = ($start>$end ? $start:$end);
	$anchors->{$anchornum}->{$accession}->{'orient'} = $orientation;
	$seq2anchors->{$accession}->{$anchornum}++;
    }
    else{
	
    }
}
#Foreach sequence, sort anchors by coordinate and print distance between adjacent coords
foreach my $accession (sort {$a cmp $b} (keys %$seq2index)){
    my @sortedanchors =  sort {$anchors->{$a}->{$accession}->{'start'} <=> $anchors->{$b}->{$accession}->{'start'}} (keys %{$seq2anchors});
    for(my $i=0;$i<scalar(@sortedanchors)-1;$i++){
	my $a1 =  $sortedanchors[$i];
	my $a2 =  $sortedanchors[$i+1];
	my $dist = &getDistance($anchors->{$a1},$anchors->{$a2});
	print STDERR "Bad coords a1:$anchors->{$a1}->{'start'} - $anchors->{$a1}->{'end'} a2:$anchors->{$a2}->{'start'} - $anchors->{$a2}->{'end'}" if($dist < 0);
	print $a1," ",$a2," ",                             #Anchors
	    $seq2index->{$accession}," ",                  #Seqindex
	    $dist," ",                                     #Distance between anchors
	$genome2index->{$accession}," ",                   #Genomeindex
	$anchors->{$a1}->{$accession}->{'orient'}," ",$anchors->{$a2}->{$accession}->{'orient'}," ", #Orientation
	$anchors->{$a1}->{$accession}->{'start'}," ",$anchors->{$a2}->{$accession}->{'start'}," ",   #Anchor1 coords
	$anchors->{$a1}->{$accession}->{'end'}," ",$anchors->{$a2}->{$accession}->{'end'}," ",       #Anchor2 coords
	"\n";
    }
}


sub getDistance{
    my($anchors1,$anchors2) = @_;
    if($anchors1->{'orient'} eq '-' && $anchors2->{'orient'} eq '-'){
	# <e----s|  <e----s|
	return $anchors2->{'end'} - $anchors1->{'start'};
    }
    elsif($anchors1->{'orient'} eq '-' && $anchors2->{'orient'} eq '+'){
	# <e---s| |s---e>
	return $anchors2->{'start'} - $anchors1->{'start'};
    }
    elsif($anchors1->{'orient'} eq '+' && $anchors2->{'orient'} eq '-'){
	# |s---e> <e---s|
	return $anchors2->{'end'} - $anchors1->{'end'};
    }
    elsif($anchors1->{'orient'} eq '+' && $anchors2->{'orient'} eq '-'){
	# |s---e> |s---e>
	return $anchors2->{'start'} - $anchors1->{'end'};
    }
    else{
	die "Bad orientations $anchors1->{'orient'} && $anchors2->{'orient'}";
	return -1;
    }

}
