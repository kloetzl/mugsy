#!/usr/bin/perl

use strict;
use AlignmentTree;
use Data::Dumper;

print STDERR "Reading $ARGV[0]\n";
my $atree = AlignmentTree::deserialize($ARGV[0]);

print STDERR "Querying $ARGV[1],$ARGV[2],$ARGV[3]\n";
my @results = $atree->intersect($ARGV[1],$ARGV[2],$ARGV[3]);

my $outputtable = [];
my $rowlookup;
my $columnlookup;

$columnlookup->{$ARGV[1]} = 1;
my $row=0;
my $column=2;
foreach my $r (@results){
    if(!exists $rowlookup->{$r->[0]}){
	$rowlookup->{$r->[0]}=$row++;
    }
    if(!exists $columnlookup->{$r->[1]}){
	$columnlookup->{$r->[1]}=$column++;
    }
}

foreach my $r (sort {$a->[2] <=> $b->[2]} @results){
    $outputtable->[$rowlookup->{$r->[0]}]->[$columnlookup->{$r->[1]}] = "$r->[2] $r->[3]";
    $outputtable->[$rowlookup->{$r->[0]}]->[0] = $r->[0];
}

my $columnwidth=20;
my $printformat='%-'.$columnwidth.'.'.$columnwidth.'s';
printf("$printformat\t","matchname");
foreach my $col (sort {$columnlookup->{$a} <=> $columnlookup->{$b}} keys %$columnlookup){
    printf("$printformat\t","$col");
}
print "\n";
foreach my $row (sort {
    if( $a->[1] eq $b->[1]){
	$b->[1] cmp $a->[1];
    }
    else{
	$a->[1] <=> $b->[1];
    }
}
		 @$outputtable){
    foreach my $col (@$row){
	$col = '-' if(!$col);
	printf("$printformat\t","$col");
    }
    print "\n";
}



