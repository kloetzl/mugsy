#!/usr/bin/perl

use strict;

my $members=0;
my $label=0;
my $blockopen=0;
my @lines=0;

while(my $line=<STDIN>){

    if($line =~ /^a score/){
	if($blockopen==1){
	    if($members>=1){
		&labelblocks(\@lines,$members,++$label);
	    }
	    else{
		print @lines;
	    }
	}
	$blockopen=1;
	$members=0;
	@lines=();
    }
    if($blockopen==1){
	push @lines,$line;
    }
    else{
	print $line;
    }
    if($line =~ /^s\s\S/){
	$members++;
    }
}
if($blockopen==1){
    if($members>1){
	&labelblocks(\@lines,$members,++$label);
    }
    else{
	print @lines;
    }
}

sub labelblocks{
    my($lines,$nummembers,$label) = @_;
    die if($lines[0] !~ /^a score=/);
    chomp $lines[0];
    my @orients;
    for(my $i=1;$i<@lines;$i++){
	my($orient) = ($lines[$i] =~ /^s\s+\S+\s+\d+\s+\d+\s+([+-])/);
	push @orients,$orient if(defined $orient);
    }
    $lines[0] .= " label=$label ";
    if($lines[0] !~ /orient/){
	$lines[0] .= " orient=+ ";
    }
    $lines[0] .= "\n";
    print @lines;
}
