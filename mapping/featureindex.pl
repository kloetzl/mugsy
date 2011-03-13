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
use lib '/usr/local/projects/angiuoli/mugsy_trunk/mapping';
use lib './';
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
    print `bp_genbank2gff3.pl --filter misc_feature -in stdin -out - < | grep -v "# Input" >> /tmp/$$.gff`;
    open FILE,"/tmp/$$.gff";
    &parseGFF(\*FILE,'gene','pseudogene');
    close FILE;
}
elsif(lc($filetype) =~ /ptt/){
    my $seqname;
    while(my $line=<STDIN>){
	if($line =~ /^>/ || $line =~ /^Location/){
	    if($line =~ /^>(\S+)/){
		$seqname = $1;
	    }
	}
	else{
	    #36..1   -       35      XOCORF_0001     -       hypothetical protein
	    my @elts = split(/\t/,$line);
	    my ($fmin,$fmax) = ($elts[0] =~ /(\d+)\.\.(\d+)/);
	    ($fmin,$fmax) = ($fmax < $fmin) ? ($fmax,$fmin) : ($fmin,$fmax);
	    $fmin = $fmin-1;
	    my $strand = $elts[1];
	    my $featname = $elts[3];
	    print "Adding feature $featname on sequence:$seqname $fmin,$fmax,$strand to alignment tree\n";
	    $atree->insert([[$seqname,$fmin,$fmax,$strand,$fmax-$fmin."M"]],'gene:'.$featname,'gene');
	    
	}
    }
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
    my $features={};
    while(my $line=<$file>){
	if($line !~ /^\#/){
	    chomp $line;
	    my @elts = split(/\t/,$line);
	    if(length($line)>0 && scalar(@elts)==9){
		if(exists $featlookup{lc($elts[2])}){
		    my %attrs = map {split(/=/)} split(/;/,$elts[8]);
		    my $geneid;
		    if(exists $attrs{'locus_tag'}){
			$geneid=$attrs{'locus_tag'};
		    }
		    elsif(exists $attrs{'ID'}){
			#Can't expect that ID is unique across files, so append sequence name
			$geneid=$elts[0].'_'.$attrs{'ID'};
		    }
		    else{
			print STDERR "Skipping unrecognized GFF3 line $line\n";
			next;
		    }
		    my $fmin = $elts[3];
		    my $fmax = $elts[4];
		    my $orient = $elts[6];
		    my $i=0;
		    
		    while(exists $features->{$geneid}){
			print "Duplicate named feature $geneid. Renaming to ${geneid}_$i\n";
			$geneid=$geneid.'_'.++$i;
		    }
		    $features->{$geneid}++;
		    die "Unsupported $fmax>=$fmin. Line: $line" if($fmax<=$fmin);
			die "Bad orient $orient. Line: $line" if($orient ne '+' && $orient ne '-');
		    $atree->insert([[$elts[0],$fmin-1,$fmax,$orient,($fmax-$fmin+1)."M"]],'gene:'.$geneid,'gene');
		}
	    }
	}
    }

}
