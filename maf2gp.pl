#!/usr/bin/perl

use strict;


print "0 0 0\n";
print "0 0 0\n";
print "\n\n";
my @x;
my $regex = ($ARGV[0] eq '+') ? '\-' : '\+';

while(my $line=<STDIN>){
    if($line =~ /^a/){
	my @p;
	foreach my $elt (@x){
	    my($acc) = ($elt =~ /s\s+(\w+)/);
	    if($acc eq $ARGV[1]){
		$p[0] = $elt;
	    }
	    if($acc eq $ARGV[2]){
		$p[1] = $elt;
	    }
	}

	&printpair(@p) if(scalar(@p) ==2);
    	@x = ();
    }
    else{
	if($line =~ /^(s.+)\s+\S+/){
	    push @x,$1;
	}
    }
}

my @p;
foreach my $elt (@x){
    my($acc) = ($elt =~ /s\s+(\w+)/);
    if($acc eq $ARGV[1] || $acc eq $ARGV[2]){
	push @p,$elt;
    }
}
&printpair(@p) if(scalar(@p) ==2);



sub printpair{
    my($ref,$qry) = @_;
    
    my($refa,$refb,$refe,$refo,$reflen) = ($ref =~ /s\s+(\w+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+)/);
    my($qrya,$qryb,$qrye,$qryo,$qrylen) = ($qry =~ /s\s+(\w+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+)/);
    $refe = $refb + $refe;
    $qrye = $qryb + $qrye;
    print "#$ref\n";
    print "#$qry\n";
    if($refo eq '+' && $qryo eq '+' && $ARGV[0] ne '-'){
	print "$refb $qryb 100\n";
	print "$refe $qrye 100\n\n\n";
    }
    elsif($refo eq '+' && $qryo eq '-' && $ARGV[0] eq '-'){
	$qrye = $qrylen - $qrye;
	$qryb = $qrylen - $qryb;
	print "$refe $qrye 100\n";
	print "$refb $qryb 100\n\n\n";

    }
    elsif($refo eq '-' && $qryo eq '+' && $ARGV[0] eq '-'){
	$refe = $reflen - $refe;
	$refb = $reflen - $refb;
	print "$refe $qrye 100\n";
	print "$refb $qryb 100\n\n\n";
    }
    else{
#		    print STDERR "$ref\n$qry\n";
    }
}
