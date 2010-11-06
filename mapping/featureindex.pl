#!/usr/bin/perl
#
#./featureindex.pl mugsyindex < mugsy.out
#Converts GFF or simple tab text files to
#Supports Genbank files if Bioperl is also installed
#
#TODO
#POD::usage
#Add more supported types from bioperl, remote download of accessions etc


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

my $filetype = $ARGV[1];

if(lc($filetype) =~ /gff/){
    print STDERR "Reading filetype $filetype\n";
    &parseGFF(\*STDIN,'gene','pseudogene');
}
elsif(lc($filetype) =~ /genbank/){
    print STDERR "Reading filetype $filetype\n";
    my $file;
    print `bp_genbank2gff3.pl --filter misc_feature -in stdin -out /tmp/$$.gff`;
    open FILE,"/tmp/$$.gff";
    &parseGFF(\*FILE,'gene','pseudogene');
    close FILE;
}
else{
    while(my $line=<STDIN>){
	my($featname,$seqname,$fmin,$fmax,$strand) = split(/\s+/,$line);
	$atree->insert([[$seqname,$fmin,$fmax,$strand,$fmax-$fmin."M"]],'gene:'.$featname,'gene');
	print "Adding feature $featname on sequence:$seqname $fmin,$fmax,$strand to alignment tree\n";
    }
}

print STDERR "Writing index to $ARGV[0]\n";
$atree->serialize($ARGV[0]);

sub parseGFF{
    my $file = shift;
    my @feattypes = @_;
    my %featlookup = map {lc($_) => 1} @feattypes;
    while(my $line=<$file>){
	if($line !~ /^\#/){
	    chomp $line;
	    my @elts = split(/\t/,$line);
	    if(length($line)>0 && scalar(@elts)==9){
		if(exists $featlookup{lc($elts[2])}){
		    my %attrs = map {split(/=/)} split(/;/,$elts[8]);
		    if(exists $attrs{'locus_tag'}){
			my $fmin = $elts[3];
			my $fmax = $elts[4];
			my $orient = $elts[6];
			
			die "Unsupported $fmax>=$fmin. Line: $line" if($fmax<=$fmin);
			die "Bad orient $orient. Line: $line" if($orient ne '+' && $orient ne '-');
			$atree->insert([[$elts[0],$fmin-1,$fmax,$orient,($fmax-$fmin+1)."M"]],'gene:'.$attrs{'locus_tag'},'gene');
		    }
		}
	    }
	}
    }

}
