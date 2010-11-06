#!/usr/bin/perl
#
#./mugsyindex.pl index.file < mugsy.out
#Adds MUGSY output to a MUGSY formatted index
#Each block is saved as type 'syntenyblk'

use strict;
use lib '/usr/local/projects/angiuoli/developer/sangiuoli/mugsy/trunk';
use AlignmentTree;
use Data::Dumper;

my $atree = new AlignmentTree();
if(-e $ARGV[0]){
    $atree = AlignmentTree::deserialize($ARGV[0]);
}
else{

}

my $currscore;
my $block = [];
my $k=0;
my $name;
while(my $line=<STDIN>){
    chomp $line;
    if($line !~ /^[\s\#]/){
	my @elts = split(/\s+/,$line);
	if($name ne $elts[0]){
	    $atree->insert($block,$name,"synteny") if(scalar @$block>0 && $name);
	    $name = "$elts[0]";
	    $block = [];
	}
	push @$block,[$elts[1],$elts[3],$elts[4],$elts[2]];
    }
}
$atree->insert($block,$name,"synteny") if(scalar @$block>0 && $name);
print STDERR "Writing index to $ARGV[0]\n";
$atree->serialize($ARGV[0]);
