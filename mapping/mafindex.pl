#!/usr/bin/perl
#
#./mafindex.pl mugsyindex < mugsy.out
#Adds an MAF formatted file to a MUGSY formatted index
#Each alignment is saved as type 'alignment'
#


use strict;
use lib '/usr/local/projects/angiuoli/developer/sangiuoli/mugsy/trunk/mapping';
use AlignmentTree;
use Storable qw(store retrieve);
use Data::Dumper;

$Storable::Deparse = 1;
$Storable::Eval = 1;

my $atree = new AlignmentTree();
if(-e $ARGV[0]){
    $atree = AlignmentTree::deserialize($ARGV[0]);
}
else{

}

my $currscore;
my $block = [];
my $k=0;
my $i=0;
my $label;
while(my $line=<STDIN>){
    if($line =~ /^a\s+score=([\d\.\-]+)/){
	my $name = "WGA_$label";
	if(exists $atree->{_alignments}->{"WGA_$label"}){
	    print "Creating new alignment name. $name taken\n";
	    $name = "WGA_".$$."_$i";
	}
	if(scalar(@$block)){
	    print "Saving alignments $name with ",scalar(@$block)," sequences\n";
	    $atree->insert($block,"$name","alignment") if(scalar(@$block));
	}
	($label) = $line =~ /label=(\w+)/;
	$label = "nolabel".++$k if !$label;
	$currscore=$1;
	$block=[];
	$i++;
    }

    elsif($line =~ /^s\s+/){
	my @elts = split(/\s+/,$line);
	#$elts[1] =~ s/\./_/g;
	#$elts[1] =~ s/\|/_/g;
	#[1] - accession
	#[2] - start
	#[3] - length
	#[4] - orient
	#[5] - seqlen
	#[6] - seq
#From UCSC FAQ about MAF format
#  start -- The start of the aligning region in the source sequence. This is a zero-based number. If the strand field is '-' then this is the start relative to the reverse-complemented source sequence.
# size -- The size of the aligning region in the source sequence. This number is equal to the number of non-dash characters in the alignment text field 
	my $start = $elts[2];
	my $end = $start+$elts[3];
	my $orient = $elts[4];
	if($orient eq '-'){
	    $start = ($elts[5] - $start - $elts[3]);
	    $end = $start + $elts[3];
	}
	my ($cigar,$len) = &get_cigar($elts[6]);
	my $seq = $elts[1];
	#Check for species.accession formatted names, trim to accession if the same
	my($species,$accession) = ($seq =~ /(\S+)\.(\S+)/);
	if($species ne "" && $species eq $accession){
	    $seq = $accession;
	}
	die "Bad orient: $orient\n" if($orient ne '-' && $orient ne '+');
	print "$seq $start,$end ", $end-$start,"\n";
	push @$block,[$seq,$start,$end,$orient,$cigar];
    }
}
my $name = "WGA_$label";
if(exists $atree->{_alignments}->{"WGA_$label"}){
    $name = "WGA_".$$."_$i";
}
print "Saving alignments $name with ",scalar(@$block)," sequences\n";
$atree->insert($block,"$name","alignment") if(scalar(@$block));
print STDERR "Writing index to $ARGV[0]\n";
$atree->serialize($ARGV[0]);


sub get_cigar{
    my($seqs) = @_;
    my $cig;
    my $len=0;
    my @chars = split(//,$seqs);
    my $count=0;
    my $curr=0; #1 - match, 2 - gap
    foreach my $c (@chars){
	#match char
	if($c ne '-'){
	    if($curr==2){
		#in gap
		#write prev gap
		$cig .= $count."X";
		$count=0;
	    }
	    #in match
	    $count++;
	    $curr=1;
	}
	else{
	    #gap char
	    if($curr==1){
		#in match
		#write prev gap
		$cig .= $count."M";
		$len += $count;
		$count=0;
	    }
	    #in gap
	    $count++;
	    $curr=2;
	}
    } 
    if($curr==1){
	#in gap
	#write prev gap
	$cig .= $count."M";
	$len += $count;
    }
    if($curr==2){
	#in gap
	#write prev gap
	$cig .= $count."X";
    }
    return ($cig,$len);
}

