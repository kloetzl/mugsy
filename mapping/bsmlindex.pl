#!/usr/bin/perl

use strict;

use XML::Twig;
use AlignmentTree;
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

my $mapping;

if(-e $ARGV[1]){
    #parsing lookup file
    open FILE, $ARGV[1] or die "Can't open mapping file $ARGV[1]";
    while(my $line=<FILE>){
	my($tseq,$oseq,$offset) = split(/\s+/,$line);
	$mapping->{$oseq} = [$tseq,$offset-1];
    }
    close FILE;
}

my $twig = new XML::Twig(
			 twig_handlers =>         
			 { 'Feature[@class = "polypeptide"]' => sub {
			     my( $twig, $elt)= @_;
			     my $iloc  = $elt->first_child('Interval-loc');
			     my $seqname  = $elt->parent('Sequence')->{'att'}->{'id'};
			     my $featname = $elt->{'att'}->{'id'};
			     my $class = $elt->{'att'}->{'class'};
			     my $complement = $iloc->{'att'}->{'complement'};
			     if ($complement eq '1'){
				 $complement = '-';
			     }
			     if ($complement eq '0'){
				 $complement = '+';
			     }
			     my($fmin,$fmax) = ($iloc->{'att'}->{'startpos'},$iloc->{'att'}->{'endpos'});
			     if(exists $mapping->{$seqname}){
				 $fmin = $fmin+$mapping->{$seqname}->[1];
				 $fmax = $fmax+$mapping->{$seqname}->[1];
				 print "Using mapping for $seqname $mapping->{$seqname}->[0]\n";
				 $seqname = $mapping->{$seqname}->[0];
			     }
			     $atree->insert([[$seqname,$fmin,$fmax,$complement,$fmax-$fmin."M"]],$featname,$class);
			     #print "$seqname\t$featname\t$fmin\t$fmax\t$class\n";
			     
			 },
		       },
			 );            

print STDERR "Writing index to $ARGV[0]\n";
$atree->serialize($ARGV[0]);

my $stdin_fh = \*STDIN;
$twig->parse($stdin_fh);
