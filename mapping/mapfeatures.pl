#!/usr/bin/perl
######################
=head1 NAME

mapfeatures - derives a set of mapped features according to a
multiple sequence alignment. Reports on the consistency of
annotated features in the mapping.

=head1 USAGE

mapfeatures.pl alignments.index seqs.fasta < features.txt 

Outputs are a series of text reports and an HTML report that can be
loaded in a web browser

Inputs:

(1) alignment.index - An index file containing a whole genome multiple
alignment and genome annotations. This index can be generated with a
combination of featureindex.pl,mafindex.pl,xmfaindex.pl. The whole
genome multiple alignment can be produced by a whole genome aligner
like Mugsy, TBA (indexed using mafindex.pl) or Mauve (index using
xmfaindex.pl). The genome annotations in Genbank or GFF3 format can be
indexed with featureindex.pl

(2) seqs.fasta - Multi-FASTA file of the input genomes. These must be
    the same genomes aligned.

(3) features.txt - A space delimited file consisting of 
feature_id sequence_id fmin fmax strand

=head1 SYNOPSIS
#############
#Example usage
#############

#Generate whole genome alignment
mugsy --prefix nmen_v16 v16/*.fsa

#Index output
mafindex.pl nmen.index < nmen_v16.maf 

#Index annotations
featureindex.pl n16.index genbank < nmen_v16.all.gbk > v16annotations.out
cat v16/*.fsa > v16.all.fsa

#Run mugsy-annotator
mugsy-annotator ./n16.index ./v16.all.fsa < v16annotations.out > v16.features.mapping

#For more detailed output (v16.html, v16.aln.report, v16.table, v16.clusters, v16.edits)
mugsy-annotator --prefix v16 --print-alignments ./n16.index ./v16.all.fsa < v16annotations.out > v16.features.mapping

#############
#APPLICATIONS
#############

1)Reporting orthologs using whole genome alignment

The script can be used to produce a list of orthologous genes in the
case where the input alignments correspond to orthologous regions

2)Reporting annotation inconsistencies, such as frameshifts or
varying start sites
Aligned annotations are further classified and checked for
consistency of start and stop codons.  Inconsistencies may indicate
annotation error, sequencing errors, or frameshifts.  Alternatively,
the inconsistencies can be due to poor or missing alignments. The
summary information provided at the end of the output provides an
indication of the overall consistency of the annotations in the set.
The script has been used to evaluate consistency of annotations
across numerous sequenced strains of bacteria and identify likely
errors

#####################
#PREPARATION OF INPUT
#####################
Meant to be used in conjunction the several utility scripts to
identify orthologs and classify annotations in a set of aligned genomes
Related scripts
>mugsymapper aln.maf features.txt > clusters.out
>mafindex
>featureindex
>gb2annottab genome.gbk1,...,genome.gbkN > orig.annot.tab
>indextab alignments.index
>updategb genome.gbk,....,genome.gbkN < annot.updates.tab

#####################
#BUGS/LIMITS/TODO
#####################

-printAlignments displays wrong frame for gene fragments that have
more than one start or end in a single display line
-Will report coverage,identity>1 if there are overlapping alignments
-Does not detect cases where gene fragments run off end of the contig
-no command line usage,help
-need to rename alignmenttree, AlignedIntervalTree

#NOTES
##########
#Input coordinates are zero start, interbase coordinates
#0 1 2 3 4
# A T A C
#The feature TA above has coordinates 1-3 
#specified in the code as fmin=1 fmax=3. Length is fmax-fmin=2
#
#Contact: S. Angiuoli (angiuoli@cs.umd.edu) 
#December 2010

=cut


use strict;
use lib '/usr/local/projects/angiuoli/mugsy_trunk/mapping';
use lib './';

use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;

#Bioperl is used only for translation machinery
use Bio::Perl;
#use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use Bio::Seq::EncodedSeq;
#use Bio::LiveSeq::Mutation; tried this but couldn't get to work properly
#Default cutoffs

use AlignmentTree;

my %options;
my $results = GetOptions (\%options, 
			  'prefix=s',
			  'input_file=s',
			  'map_file=s',
			  'featlist=s', #Restrict mapping to list of features
			  'duplications=s', #Report duplications, requires addl index file
			  'coverage|c=s',
			  'query_coverage|q=s',
			  'identity|i=s',
		
			  'sortkeys=s', #sort order of fields 'gfreq','len','afreq' when reporting edits
			  'reportedits=s', #top number of edits to report, default all
			  'maxchange=s', #max allowable %length changes
			  'prefix=s', #Generate output reports with file prefix
			  'cogformat=s', #Output cog format to stdout
			  'printalignments', 
			  'printhtml',
			  'skipframeshifts',
			  #Missing gene options
			  'minorflen=s',
			  'maxorflen=s',
			  'verbose|v', #Verbose warnings
			  'debug|d=s') || pod2usage(-verbose => 1);

pod2usage(-verbose=>1) if($options{'help'});


my $coverage_cutoff = (exists $options{'coverage'}) ?  $options{'coverage'} : 0.5;
my $query_coverage_cutoff = (exists $options{'query_coverage'}) ?  $options{'query_coverage'} : 0;
my $pid_cutoff= (exists $options{'identity'}) ?  $options{'identity'} : 0.1;
print "#Using coverage cutoff:$coverage_cutoff identity:$pid_cutoff query_coverage:$query_coverage_cutoff\n";

my $MAXORFLEN = (exists $options{'maxorflen'}) ? $options{'maxorflen'} : 30000; #in bp
my $MINORFLEN = $options{'minorflen'} || 30; #in aa residues
my $ORFLEN_MAXDELTA = 0.5; #do not consider possible codons that are less than X the length of the maximum annotated ORF

my $FS_THRESHOLD = 3;
my $FSLEN_THRESHOLD = 9;
#'frameshift_consistency=s', #only report frameshifts that occur in < X fraction of aligned sequences. 1 show all possible frameshifts, default 0.5.
#my $FS_FRACTIONGENOME = .5;
#my $FS_FRACTIONGENOME = (exists $options{'frameshift_consistency'}) ? $options{'frameshift_consistency'} : 0.5;;



#Used for detecting contig boundaries
my $PMARK_SPACER = "NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN";

#Flag for checking consistent start,stop
#Assumes input features are genes
my $doconsistencychecks=1;

#Report new ORFs using aligned start codons
my $dofindneworfs = 0;
my $autocorrect=0;

#Only report alternative start codons that
#results in a longer ORF
my $longer_altstarts=1;
my $moreconsistent_altstarts=1;

#Only report alternative start codons that
#appear more frequently in the aligned genoems
my $freq_altstarts=1;
my $freq_altstops=0;

my $aligntoken="WGA";
my $CODON_DELIM = '.';
my $CODON_DELIM_REGEX = '\.';

#Output flags
my $COGoutputformat=(exists $options{'cogformat'}) ? $options{'cogformat'} : 0;
my $cogfh; #cog format
my $cfh; #cluster format
my $ctfh; #table
my $ctfh2; #table with coords
my $htmlout=(exists $options{'printhtml'} ? 1 : 0);

if($COGoutputformat){
    open $cogfh, "+>$options{'cogformat'}" or die "Can't open COG file $options{'cogformat'}";#\*STDOUT;
}
else{
    open $cogfh, "+>$options{'prefix'}clusters.cog";
}
if(! $options{'prefix'}){
    $options{'prefix'} = "mugsyant/run$$";
    print `mkdir -p $options{'prefix'}`;
}
print STDERR "Writing output to $options{'prefix'}clusters.table, $options{'prefix'}clusters.coords.table, $options{'prefix'}clusters.out \n";
print "#Writing output to $options{'prefix'}clusters.table, $options{'prefix'}clusters.coords.table, $options{'prefix'}clusters.out \n";
print "#EDITTBL format cluster_id, codon_id, genome_freq, currentannotated_freq, avglen, num_orgswithoverlaps, comments\n";
open $cfh, "+>$options{'prefix'}clusters.out";
open $ctfh, "+>$options{'prefix'}clusters.table";
open $ctfh2, "+>$options{'prefix'}clusters.coords.table";


my $printalignments=(exists $options{'printalignments'}) ? $options{'printalignments'} : 0;
my @sortkeys = (exists $options{'sortkeys'}) ? (split(/,/,$options{'sortkeys'})) : ('gfreq','len','afreq');
if(scalar @sortkeys != 3){

    print STDERR "Enter sort order using names gfreq,afreq,len for aligned frequency in the genome, annotated frequency, and ORF length. Sort is in descending order, largest value first.\n";
    print STDERR "eg. --sortkeys gfreq,len,afreq\n";

    exit 1;
}

#Debugging flags
my $checkbadlen=0;
my $debug=$options{'debug'};
my $printskipped=$debug;
my $verbose=(exists $options{'verbose'} ? 1 : 0); #verbose warnings, mostly for debugging
if($debug){
    $verbose=$debug;
}

#Master list of features and attributes
#0-seqname
#1-fmin
#2-fmax
#3-len
#4-orient
#5-polyid
#6-geneid
#7-startcodon pos
#8-startcodon aln
#9-stopcodon pos
#10-stopcodon aln
my $features = {};
my $allseqs = {};
my $seqindex = {};

#Featlist
my $featlist;
if(exists $options{'featlist'}){
    foreach my $f (split(/\s+/,$options{'featlist'})){
	$featlist->{$f}=1;
    }
}
if($featlist && scalar (keys %$featlist)){
	print STDERR "Limiting results to ",scalar(keys %$featlist)," genes\n";
}
my $codons = {};
my $classes_sum = {};
my $classes_all = {};

my $newclasses_sum = {};

#AlignmentTree is a interval tree that contains alignments between sequences
#and features on those sequences
my $atree = AlignmentTree::deserialize($ARGV[0]);
$atree->{_debug}=$debug;
#Read a white space delimited list of features to map from stdin

my $datree;
if(-e $options{'duplications'}){
    $datree = AlignmentTree::deserialize($options{'duplications'});
}

my %featlookup;
my $filetype;
my $fh;

$options{'input_file'} = $options{'map_file'} if(exists $options{'map_file'});

if($options{'input_file'}) {
	open($fh, "<$options{'input_file'}") or die "Error in opening the file, $options{'input_file'}, $!\n";
} else {
	$fh = \*STDIN;
}
my $seqref;
while(my $line=<$fh>){
    my($name,$seq,$fmin,$fmax,$orient,$polyid,$geneid,$annotations);
    chomp $line;
    if($line =~ /\#gff-version 3/){
	$filetype = 'gff3';
	$featlookup{'gene'}++;
	$featlookup{'pseudogene'}++;
    }
    elsif($line =~ /^>/ || $line =~ /^Location/){
	$filetype = 'ptt';
	if($line =~ /^>(\S+)/){
	    $seqref = $1;

	}
    }
    elsif($line !~ /^\#/){
	if($filetype eq 'gff3'){
	    #GFF
	    my @elts = split(/\t/,$line);
	    if(scalar(@elts)==9){
		if(exists $featlookup{lc($elts[2])}){
		    my %attrs = map {split(/=/)} split(/;/,$elts[8]);
		    if(exists $attrs{'locus_tag'}){
			$name = $attrs{'locus_tag'};
		    }
		    elsif(exists $attrs{'ID'}){
                        #Can't expect that ID is unique across files, so append sequence name
                        $name=$elts[0].'_'.$attrs{'ID'};
                    }

		    if(exists $attrs{'product'}){
			$annotations .= $attrs{'product'};
		    }
		    if(lc($elts[2]) eq 'pseudogene'){
			$annotations .= "pseudogene ";
		    }
		    ($seq,$fmin,$fmax,$orient,$polyid,$geneid) = ($elts[0],$elts[3],$elts[4],$elts[6],$name,$name);
		    ($fmin,$fmax) = ($fmin<$fmax) ? ($fmin-1,$fmax) : ($fmax-1,$fmin);
		}
		elsif(lc($elts[2]) eq 'cds'){
		    #hack for names from genbank
		    my %attrs = map {split(/=/)} split(/;/,$elts[8]);
		    if(exists $attrs{'product'}){
			$annotations .= $attrs{'product'};
			my $cdsname;
			if(exists $attrs{'locus_tag'}){
			    $cdsname = $attrs{'locus_tag'};
			    if(exists $features->{$cdsname}){
				$features->{$cdsname}->[11] = $annotations;
			    }
			}
			if(exists $attrs{'ID'} && (length($cdsname)>0 && !exists $features->{$cdsname})){
			    #Can't expect that ID is unique across files, so append sequence name 
			    $cdsname=$elts[0].'_'.$attrs{'ID'};
			    if(exists $features->{$cdsname}){
				$features->{$cdsname}->[11] = $annotations;
			    }
			}
		    }
		}
	    }
	    else{
		#print "Skipping $line\n" if($debug);
	    }
	}
	elsif($filetype eq 'ptt'){
	    #36..1   -       35      XOCORF_0001     -       hypothetical protein
	    my @elts = split(/\t/,$line);
	    ($fmin,$fmax) = ($elts[0] =~ /(\d+)\.\.(\d+)/);
	    ($fmin,$fmax) = ($fmax < $fmin) ? ($fmax,$fmin) : ($fmin,$fmax);
	    $orient = $elts[1];
	    $fmin = $fmin-1;
		
	    $name = $elts[3];
	    $seq = $seqref;
	    $annotations .= $elts[5];
	}
	else{
	    #Custom simple space delim text
	    my @annots;
	    ($name,$seq,$fmin,$fmax,$orient,$polyid,$geneid,@annots) = split(/\s+/,$line);
	    $annotations .= join (' ',@annots);
	    #Allow for 0,1 orient
	    if($orient =~ /\d/){
		if($orient > 0){
		    $orient = '+';
		}
		else{
		    $orient = '-';
		}
	    } 
	    die "Bad orient $orient\n" if($orient ne '-' && $orient ne '+');
	}
	if(length($name)>0){
	    die "Unsupported $fmax>=$fmin. $line" if($fmax<=$fmin);
	    if($fmin<0){
		print STDERR "Illegal fmin $fmin for $seq,$fmin,$fmax,$fmax-$fmin,$orient,$polyid,$geneid\n";
		$fmin=0;
		next;
	    }
	    die "Bad orient $orient. $line" if($orient ne '+' && $orient ne '-');
	    my $i=0;
	    while(exists $features->{$name}){
		print "#Duplicate named feature $name. Renaming to ${name}_$i\n";
		$name=$name.'_'.++$i;
	    }
	    if(!defined $featlist || exists $featlist->{$name}){
		$features->{$name} = [$seq,$fmin,$fmax,$fmax-$fmin,$orient,$polyid,$geneid];
		
		my($org) = ($seq =~ /([^\.]+)/);
		$allseqs->{$org}++;
		#[7]-[10] reserved for start,stop codon info
		$features->{$name}->[11] = $annotations;
	    }
	}
    }
}

my @sortedallseqs = sort {$a cmp $b} (keys %$allseqs);
for(my $i=0;$i<scalar(@sortedallseqs);$i++){
    $seqindex->{$sortedallseqs[$i]}=$i;
}

#Save a list of clusters
my $clusters = {};
#Current cluster id, a unique identifier for a cluster
my $cluster_id = 0;
#Count of clusters that pass cutoffs
my $validcluster = 0;

#All genes are categorized into one of three categories
#mapped   - aligned to other genes in the set above cutoffs
#unmapped - aligned to other genes in the set but none above cutoffs
#nomatches- not aligned to any other genes in the input set
#List of mapped,unmapped,nohit genes
my $mapped = {};
my $unmapped = {};
my $deleted = {};
my $nomatches = {};
my $dups = {};
my $neworfcount = 0;
my $adjustedorfs = 0;

#List of newly called ORFs
my $neworfs = {};
#and the annotated ORFs they replace
my $subsumed = {};

#Map of feature => organism
my $feat2organism = {};

my $db;


if(-f "$ARGV[1]"){
    print STDERR "Using FASTA file $ARGV[1]. Debug level: $debug\n";
    #Faster to read everything into RAM
    #$db = Bio::DB::Fasta->new($ARGV[1],'-reindex'=>1); 
    my @ids;# = $db->ids();
    my $istream = Bio::SeqIO->new(-file => $ARGV[1],
				  -format => 'Fasta');
    while ( my $seq = $istream->next_seq()){
	push @ids, $seq->id();
	$db->{$seq->id()} = $seq;
	print "#Storing ",$seq->id(),"\n" if($verbose);
    }
    print "#Parsed FASTA sequences for ",join(',',@ids),"\n";
}
else{
    print STDERR "No FASTA file provided. Reporting alternative start codons but not calling ORFs\n";
}


#The mapping algorithm builds clusters of aligned genes in a greedy
#fashion, starting with the longest feature in the input set and
#mapping all aligned features that pass cutoffs. In the case of where
#features are genes and the alignments are orthologous regions, such
#as those identified by whole genome alignments(WGA), the clusters
#represent orthologous genes.

#Sort query genes by length in decreasing order, longest to
#shortest. In doing do, all aligned genes that cover the query gene
#above cutoffs are considered putative orthologs to the query. And the
#query gene is always the longest member of the cluster. The reported
#%id and %cov are relative to the query

foreach my $query (sort {$features->{$b}->[3] <=> $features->{$a}->[3]} #Sort on length, decreasing order
		   keys %$features){                                    #Over all features
    print "#Processing $query ",`date` if($verbose);
    #As the algorithm progresses, features are mapped and removed from consideration
    #Consider genes that remain unmapped or 
    #remain covered by <= cutoff% of length in alignments already considered
    if(!exists $mapped->{$query} && !exists $deleted->{$query}){

	#Start a new cluster based on the query gene.  Set a new
	#cluster id; each cluster can also be identified by the query
	#gene ($query)
	$cluster_id++;
	
	my($mappedorgs,$mappedgenes,$unmappedorgs,$unmappedgenes) = &buildCluster($atree,$query);
	print "#MAPPED Num_orgs:",scalar(keys %$mappedorgs)," Num_genes:",scalar(keys %$mappedgenes)," UNMAPPED Num_orgs:",scalar(keys %$unmappedorgs)," Num_genes:",scalar(keys %$unmappedgenes),"\n" if($verbose);
	die "Less than 2 mapped sequences" if(scalar(keys %$mappedgenes)>1 && scalar(keys %$mappedorgs)<=0);
	die "No mapped genes" if(scalar(keys%$mappedgenes)<1);

	#Mark inconsistencies in the cluster and save start,stop codon
	#positions of annotated genes only
	#Codon aligned, annotated frequency is also saved as 
	#'start','stop',=>seqname
	#'pairs'
	#=>
	# 'gfreq' -aligned genomic freq 
	# 'afreq' -annotated freq
	# 'len' - average length 
	#
	my($feat_attrs,$cluster_attrs,$codons) = &annotateCluster($atree,$mappedgenes,$mappedorgs);

	my $seq_attrs = {};
	my $new_orfs = {};
	
	#Look for unannotated ORFs in remaining aligned seqs using other annotated start codons
	#This can also recall orfs in the unmapped set
	#if($dofindneworfs && !$options{'skipneworfs'}){
	#$new_orfs = &findnewORFs($db,$atree,$mappedorgs,$mappedgenes,$codons);	    
	#}
	
	if((scalar(keys %$mappedgenes)>1 && scalar(keys %$mappedorgs)>1)){
	    print "#Cluster WGA$cluster_id codon_pairs:",scalar(keys %{$codons->{'pairs'}}),"\n" if($verbose);
	
	    #We have a good cluster, save it
	    #Save the cov,pid in master list of mapped genes

	    my $totallen=0;
	    my $maxlen=0;
	    foreach my $feat_name (keys %$mappedgenes){
		die "Feature $feat_name already mapped" if(exists $mapped->{$feat_name});
		$mapped->{$feat_name}->{'cov'}=$mappedgenes->{$feat_name}->{'cov'}/$features->{$feat_name}->[3];
		$mapped->{$feat_name}->{'pid'}=$mappedgenes->{$feat_name}->{'pid'}/$mappedgenes->{$feat_name}->{'len'};
		$mapped->{$feat_name}->{'cluster_id'}=$cluster_id;
		$totallen += $features->{$feat_name}->[3];
		$maxlen = ($features->{$feat_name}->[3] > $maxlen) ? $features->{$feat_name}->[3] : $maxlen;
		delete $unmapped->{$feat_name};
		die if(exists $unmappedgenes->{$feat_name});
	    }
	    my $avglen=$totallen/(scalar keys %$mappedgenes);
	    my $classesstr;
	    my $classesallstr;

	    #Save alternative ORFs
	    my $altcodons = {};

	    if(!defined $options{'reportedits'} || $options{'reportedits'} > 0){
		#Save aligned and annotated codon frequency
		foreach my $p (keys %{$codons->{'pairs'}}){
		    print "#Analyzing codon pair $p\n" if($verbose);
		    my($startcodon,$stopcodon) = split(/:/,$p);
		    foreach my $seqname (keys %$mappedorgs,keys %$unmappedorgs){
			print "#Sequence $seqname\n" if($verbose);
			#if this is the annotated pair
			if(exists $codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} && $codons->{'pairs'}->{$p}->{'orgs'}->{$seqname}->[3]==1){
			    #Do nothing, already annotated
			    $codons->{'pairs'}->{$p}->{'features'}->{$seqname} = $mappedorgs->{$seqname}->{'features'};
			    my($fmin,$fmax,$orient) = &findCoords($atree,$seqname,$startcodon,$stopcodon);
			    my $isorf = &isORF($db,$seqname,$fmin,$fmax,$orient);
			    if($isorf<=0){
				print "#BAD ORF $seqname,$fmin,$fmax ",join(',',keys %{$mappedorgs->{$seqname}->{'features'}}),"\n" if($verbose);
				foreach my $feat (keys %{$mappedorgs->{$seqname}->{'features'}}){
				    $feat_attrs->{$feat}->{'CX'}++;
				}
			    }
			    print "#annotated\n" if($verbose);
			}
			else{
			    print "#checking\n" if($verbose);
			    #check if this is an ORF in $seqname
			    my($fmin,$fmax,$orient) = &findCoords($atree,$seqname,$startcodon,$stopcodon);
			    if(defined $fmin && defined $fmax && defined $orient && $fmax>$fmin){
				die "$atree,$seqname,$startcodon,$stopcodon" if(! defined $fmin || ! defined $fmax);
				my $isorf = &isORF($db,$seqname,$fmin,$fmax,$orient);
				if($isorf>0){
				    print "#isORF true\n" if($verbose);
				    #There is an ORF on $seqname over this interval
				    if(exists $unmappedorgs->{$seqname}){
					die if(exists $mappedorgs->{$seqname});
					#genome is aligned but no ORFs above cutoffs
					if(exists $unmappedorgs->{$seqname}->{'features'}){
					    #the region is annotated with an ORF that matches below cutoffs
					    #requires a new ORF that is different than currently annotated
					    $codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,-1];
					    $codons->{'pairs'}->{$p}->{'features'}->{$seqname} = $unmappedorgs->{$seqname}->{'features'};	
					}
					else{
					    #the region is not annotated
					    #Requires a new ORF in an unannotated region
					    my $olapgenes = &getFeaturesByInterval($atree,$seqname,$fmin,$fmax,$orient);
					    my $nummapped=0;
					    if(scalar(keys %$olapgenes)>0){
						print "#Found ",scalar(keys %$olapgenes)," in region $seqname,$fmin,$fmax with no mapped,unmapped\n" if($debug);
						foreach my $gene (keys %$olapgenes){
						    if(exists $mapped->{$gene}){
							$nummapped++;
						    }
						}
					    }
					    if(scalar (keys %$olapgenes)==0){
                                                #the region is not annotated
						#Requires a new ORF in an unannotated region
						$codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,-2];
						print "#neworf $p $seqname ",$fmax-$fmin,"\n" if($verbose);
					    }
					    elsif($nummapped==0){
						#the region is annotated
						#Requires a alt ORF in an unannotated region
						$codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,-1];
						$codons->{'pairs'}->{$p}->{'features'}->{$seqname} = $unmappedorgs->{$seqname}->{'features'};
					    }
					}
				    }
				    else{
					#the region is aligned with an annotated ORF above cutoffs 
					#requires a new ORF that is different than currently annotated
					die if(! exists $mappedorgs->{$seqname});
					$codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,0];
					if(exists $mappedorgs->{$seqname}){
					    $codons->{'pairs'}->{$p}->{'features'}->{$seqname} = $mappedorgs->{$seqname}->{'features'};	
					    print "#altorf\n" if($verbose);
					}
					else{
					    print "#altorf, prev did not pass cutoffs\n" if($verbose);
					}
				    }
				}
				elsif($isorf==0){
				    print "#isORF false\n" if($verbose);
				    if(! defined $options{'skipframeshifts'}){
					if(($fmax-$fmin)<$MAXORFLEN && ($fmax-$fmin)>$MINORFLEN){
					    #See if we can call an ORF over this region with frameshifts
					    my $feat_name;
					    my $annotatedstop;
					    my $annotatedstart;
					    if(exists $mappedorgs->{$seqname} && scalar(keys %{$mappedorgs->{$seqname}->{'features'}}) == 1){
						$feat_name = [keys %{$mappedorgs->{$seqname}->{'features'}}]->[0];
						$annotatedstop = $features->{$feat_name}->[9] . $CODON_DELIM . $features->{$feat_name}->[10];
						$annotatedstart = $features->{$feat_name}->[7] . $CODON_DELIM . $features->{$feat_name}->[8];
					    }
					    #Look for possible frameshifts if there are either
					    if(exists $unmappedorgs->{$seqname}                             #a) No annotations aligned in this region above cutoffs
					       || 
					       (exists $mappedorgs->{$seqname}
						&&
						(scalar(keys %{$mappedorgs->{$seqname}->{'features'}}) > 1 	#b) Multiple annotated ORFs in this region
						 
						 ||
						 $stopcodon ne $annotatedstop 			        #c) Single annotated ORF with a different stop codon
						 ||
						 $startcodon ne $annotatedstart
						 )
						)){
						die "$seqname found in both mapped and unmapped org lists" if(exists $unmappedorgs->{$seqname} && exists $mappedorgs->{$seqname});
						print "#Considering FS on $seqname for pair $startcodon,$stopcodon annotated:$annotatedstart,$annotatedstop $fmin,$fmax,$orient\n" if($verbose);
						print "#Considering FS $stopcodon ne $annotatedstop for $feat_name on $seqname\n" if($debug && exists $mappedorgs->{$seqname});
						#Find most similar sequence that has this ORF
						my($nearestseq) = &findNearestNeighbor($atree,$seqname,$mappedorgs,$fmin,$fmax);
						print "#Using $nearestseq as nearest neighbor to $seqname\n" if($verbose);
						#Look for frameshifting mutations in $seqname
						my($fs,$netfs) = &reportFrameShifts($atree,$db,$seqname,$nearestseq,$startcodon,$stopcodon);
						if(ref $fs){
						    my $isorf = &isORF($db,$seqname,$fmin,$fmax,$orient,$fs);
						    print "#Possible ORF with frameshift indels:",scalar(@$fs)," net:$netfs isorf:$isorf\n" if($verbose);
						    if($isorf>0){
							if($isorf==2 || abs($netfs) < $FS_THRESHOLD){
							    foreach my $fs (sort {$a->[0] <=> $b->[0]} @$fs){
								if($verbose){
								    print "#FS ",join(',',@$fs)," $netfs $isorf\n";
								}
							    }
							    print "#Adding frameshift net:",scalar(@$fs)," $netfs\n" if($debug);
							    if(exists $mappedorgs->{$seqname}){
								$codons->{'pairs'}->{$p}->{'features'}->{$seqname} = $mappedorgs->{$seqname}->{'features'};
								$codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,0,$fs];
							    }
							    else{
								my $olapgenes = &getFeaturesByInterval($atree,$seqname,$fmin,$fmax,$orient);
								
								if(scalar(keys %$olapgenes)>0){
								    $codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,0,$fs];
								}
								else{
								    print "#Neworf in frameshifted region $seqname,$fmin,$fmax,$orient\n" if($debug);
								    $codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,-2,$fs];
								}
							    }
							}
						    }
						}
					    }
					}
					else{
					    print "#Skipping frameshifts check on $seqname range too big or small $fmin-$fmax\n" if($verbose);
					}
				    }
				}
				elsif($isorf == -1){
				    #TODO, check point mutation of start, stop codon
				}
			    }
			}
		    }
		}
		foreach my $p (keys %{$codons->{'pairs'}}){
		    my $fscount=0;
		    my $orgcount = scalar(keys %{$codons->{'pairs'}->{$p}->{'orgs'}});
		    foreach my $org (keys %{$codons->{'pairs'}->{$p}->{'orgs'}}){
			if(ref $codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[4]){
			    $fscount++;
			}
		    }
		    foreach my $org (keys %{$codons->{'pairs'}->{$p}->{'orgs'}}){
			#Skip frameshift if more occurs in more than $FS_FRACTIONGENOME
			if(ref $codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[4]){
			    if(1){#$fscount/$orgcount <= $FS_FRACTIONGENOME){
				$codons->{'pairs'}->{$p}->{'fsvars'}->{$org}=$codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[4];
				print "#FS ",join(',',@{$codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[4]}),"\n" if($debug);
			    }
			    else{
				next;
			    }
			}
			print "#CODONPAIR ",join(',',@{$codons->{'pairs'}->{$p}->{'orgs'}->{$org}}),"\n" if($debug);
			
			$codons->{'pairs'}->{$p}->{'gfreq'}++;
			$codons->{'pairs'}->{$p}->{'afreq'}++ if($codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[3]==1); #inc only if annotated
			$codons->{'pairs'}->{$p}->{'length'}+=($codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[1] - $codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[0]);
			#-2 encodes a new orf, no prior annotation on this $org
			if($codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[3] == -2){
			    $codons->{'pairs'}->{$p}->{'neworfs'}->{$org}->{'fmin'} = $codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[0];
			    $codons->{'pairs'}->{$p}->{'neworfs'}->{$org}->{'fmax'} = $codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[1];
			    $codons->{'pairs'}->{$p}->{'neworfs'}->{$org}->{'orient'} = $codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[2];
			}
		    }
		    $codons->{'pairs'}->{$p}->{'len'} = $codons->{'pairs'}->{$p}->{'length'}/$codons->{'pairs'}->{$p}->{'gfreq'} if($codons->{'pairs'}->{$p}->{'gfreq'} > 0);
		}
		
		$classesstr = join(';',sort {$a cmp $b} keys %{$cluster_attrs});
		#Suggest edits for inconsistently annotated clusters
		if($doconsistencychecks){
		    #Choose N best start,stop pairs according to sortkeys
		    my @bestcodonpair = sort {
			if($codons->{'pairs'}->{$a}->{$sortkeys[0]} eq $codons->{'pairs'}->{$b}->{$sortkeys[0]}){
			    if($codons->{'pairs'}->{$a}->{$sortkeys[1]} eq $codons->{'pairs'}->{$b}->{$sortkeys[1]}){
				#sort on tertiary sortkey, eg length
				$codons->{'pairs'}->{$b}->{$sortkeys[2]} <=> $codons->{'pairs'}->{$a}->{$sortkeys[2]};
			    }
			    else{
				#sort on secondary sortkey, eg annotated frequency
				$codons->{'pairs'}->{$b}->{$sortkeys[1]} <=> $codons->{'pairs'}->{$a}->{$sortkeys[1]};
			    }
			}
			else{
			    #sort on primary sortkey, eg. aligned frequency of start codon in the genome
			    $codons->{'pairs'}->{$b}->{'gfreq'} <=> $codons->{'pairs'}->{$a}->{'gfreq'};
			}
		    } (keys %{$codons->{'pairs'}});
		    if(scalar(@bestcodonpair)>0){
			open EFILE, ">$options{'prefix'}cluster$cluster_id.edits.out";
			for(my $i=0;$i<scalar(@bestcodonpair);$i++){
			    my $bestcodon = $bestcodonpair[$i];
			    if($codons->{'pairs'}->{$bestcodon}->{'gfreq'} > 1){
				my $codonlength = ($codons->{'pairs'}->{$bestcodon}->{'length'}/$codons->{'pairs'}->{$bestcodon}->{'gfreq'});
				
				my $deltafracmax = ($maxlen-$codonlength)/$maxlen;
				if($deltafracmax < $ORFLEN_MAXDELTA && $deltafracmax > (-1)*$ORFLEN_MAXDELTA){
				    
				    print EFILE ">CLUSTER_$cluster_id $bestcodon\n";
				    
				    foreach my $org (keys %{$codons->{'pairs'}->{$bestcodon}->{'orgs'}}){
					#Check if there are existing annotations
					if(scalar(keys %{$codons->{'pairs'}->{$bestcodon}->{'features'}->{$org}})>0){
					    foreach my $feat_name (keys %{$codons->{'pairs'}->{$bestcodon}->{'features'}->{$org}}){
						my $pred_feat = $codons->{'pairs'}->{$bestcodon}->{'orgs'}->{$org};
						#Check if codon pair results in an alternative ORF
						if($pred_feat->[0] ne $features->{$feat_name}->[1] || 
						   $pred_feat->[1] ne $features->{$feat_name}->[2]){
						    my $fs = (defined $codons->{'pairs'}->{$bestcodon}->{'orgs'}->{$org}->[4]) ? "F" : "";
						    #if($fs eq 'F'){
							#die if(! exists $codons->{'pairs'}->{$bestcodon}->{'fsvars'}->{$org});
							#my @fsruns = @{$codons->{'pairs'}->{$bestcodon}->{'fsvars'}->{$org}};
							#foreach my $r (@fsruns){
							#    $fs .= print "[$r->[0]-$r->[1] $r->[2]:$r->[3]] $r->[4]";
							#}
						    #}
						    my $olapgenes = &getFeaturesByInterval($atree,$org,$pred_feat->[0],$pred_feat->[1],$pred_feat->[2]);
						    my $olaps;
						    foreach my $feat (keys %$olapgenes){
							if($feat ne $feat_name && ! exists $codons->{'pairs'}->{$bestcodon}->{'features'}->{$org}->{$feat}){
							    print "#Overlapping gene found on $org ",join(",",@{$olapgenes->{$feat}}),"\n" if($verbose);
							    $olaps->{$feat} = $olapgenes->{$feat};
							    $codons->{'pairs'}->{$bestcodon}->{'olaps'}->{$org}++
							}
						    }
						    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'fmin'}=$pred_feat->[0];
						    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'fmax'}=$pred_feat->[1];
						    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'orient'}=$pred_feat->[2];
						    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'fs'}=$fs;
						    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'olaps'}=$olaps;
						    my $sdist;
						    my $sameframe=1;
						    if($pred_feat->[2] eq $features->{$feat_name}->[4]){
							$sameframe=1;
						    }
						    else{
							$sameframe=0;
						    }
						    if($pred_feat->[2] eq '+'){
							$sdist = $pred_feat->[0]-$features->{$feat_name}->[1];
						    }
						    elsif($pred_feat->[2] eq '-'){
							$sdist = $pred_feat->[1]-$features->{$feat_name}->[2];
						    }
						    print EFILE "$feat_name\t$org\t$pred_feat->[0]\t$pred_feat->[1]\t",($pred_feat->[1] - $pred_feat->[0]),"\t$pred_feat->[2]\t$sameframe\t$sdist\t$fs\t";
						    my @olaplist = keys %$olaps;
						    for(my $i=0;$i<scalar(@olaplist);$i++){
							print EFILE "$olaplist[$i]($olaps->{$olaplist[$i]}->[3];$olaps->{$olaplist[$i]}->[6] bp;";
							die "Bad feat $olaplist[$i]" if(! exists $features->{$olaplist[$i]} && scalar(keys %$featlist)==0);
							printf EFILE "%.1f)",($olaps->{$olaplist[$i]}->[6]/$features->{$olaplist[$i]}->[3]) if(exists $features->{$olaplist[$i]});

							print EFILE "," if($i<scalar(@olaplist)-1);
						    }
						    print EFILE "\n";
						}
					    }
					}
					else{
					    #Report an annotation in a region with annotations below cutoffs or neworfs
					    my $pred_feat = $codons->{'pairs'}->{$bestcodon}->{'orgs'}->{$org};
					    my $fs = (defined $codons->{'pairs'}->{$bestcodon}->{'orgs'}->{$org}->[4]) ? "F" : "";
					    my $olapgenes = &getFeaturesByInterval($atree,$org,$pred_feat->[0],$pred_feat->[1],$pred_feat->[2]);
					    
					    my $feat_name;
					    if(scalar(keys %$olapgenes)>0){
						$feat_name = "ALTORF_C$cluster_id";
					    }
					    else{
						#NEWORF due to frame
						if(! exists $codons->{'pairs'}->{$bestcodon}->{'neworfs'}->{$org}){
						    print STDERR "Unexpected neworf NEWORF_C$cluster_id: ",join(',',@$pred_feat),"\n";
						    $codons->{'pairs'}->{$bestcodon}->{'neworfs'}->{$org}->{'fmin'} = $pred_feat->[0];
						    $codons->{'pairs'}->{$bestcodon}->{'neworfs'}->{$org}->{'fmax'} = $pred_feat->[1];
						    $codons->{'pairs'}->{$bestcodon}->{'neworfs'}->{$org}->{'orient'} = $pred_feat->[2];
						}
						$feat_name = "NEWORF_C$cluster_id";
					    }
					    my $sdist;
					    my $sameframe=1;
					    if(exists $features->{[keys %$olapgenes]->[0]}){
						if($pred_feat->[2] eq $features->{[keys %$olapgenes]->[0]}->[4]){
						    $sameframe=1;
						}
						else{
						    $sameframe=0;
						}
						if($pred_feat->[2] eq '+'){
						    $sdist = $pred_feat->[0]-$features->{[keys %$olapgenes]->[0]}->[1];
						}
						elsif($pred_feat->[2] eq '-'){
						    $sdist = $pred_feat->[1]-$features->{[keys %$olapgenes]->[0]}->[2];
						}
					    }
					    #name,org,fmin,fmax,len,orient,sameframe,startdist,fs,overlaps
					    print EFILE "$feat_name\t$org\t$pred_feat->[0]\t$pred_feat->[1]\t",($pred_feat->[1] - $pred_feat->[0]),"\t$pred_feat->[2]\t$sameframe\t$sdist\t$fs\t";
					    
					    my @olaplist = keys %$olapgenes;
					    for(my $i=0;$i<scalar(@olaplist);$i++){
						print EFILE "$olaplist[$i]($olapgenes->{$olaplist[$i]}->[3];$olapgenes->{$olaplist[$i]}->[6] bp;";
						die if(! exists $features->{$olaplist[$i]} && scalar(keys %$featlist)==0);
						printf EFILE "%.1f)",($olapgenes->{$olaplist[$i]}->[6]/$features->{$olaplist[$i]}->[3]) if(exists $features->{$olaplist[$i]});
						print EFILE "," if($i<scalar(@olaplist)-1);
						$codons->{'pairs'}->{$bestcodon}->{'olapgenes'}->{$org}++;
					    }
					    print EFILE "\n";
					    if(exists $codons->{'pairs'}->{$bestcodon}->{'neworfs'}->{$org} && scalar(keys %$olapgenes)>0){
						foreach my $gene (keys %$olapgenes){
						    if(exists $mapped->{$gene}){
							print STDERR "#Neworfs marked in region on $org with genes already mapped into clusters ",join(',',@{$olapgenes->{$gene}}),"\n";
						    }
						    else{
							print STDERR "#Neworfs marked in region on $org with other annotations ",join(',',@{$olapgenes->{$gene}}),"\n";
						    }
						}
					    }
					    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'fmin'}=$pred_feat->[0];
					    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'fmax'}=$pred_feat->[1];
					    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'orient'}=$pred_feat->[2];
					    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'fs'}=$fs;
					    $altcodons->{$bestcodon}->{'orgs'}->{$org}->{'olaps'}=$olapgenes;

					}
				    }
				    

				    $altcodons->{$bestcodon}->{'name'}="ALT$i";
				    $altcodons->{$bestcodon}->{'gfreq'}=$codons->{'pairs'}->{$bestcodon}->{'gfreq'};
				    $altcodons->{$bestcodon}->{'afreq'}=$codons->{'pairs'}->{$bestcodon}->{'afreq'};
				    $altcodons->{$bestcodon}->{'len'}=$codons->{'pairs'}->{$bestcodon}->{'len'};
				    $altcodons->{$bestcodon}->{'neworfs'}=$codons->{'pairs'}->{$bestcodon}->{'neworfs'};
				    $altcodons->{$bestcodon}->{'fs'}=$codons->{'pairs'}->{$bestcodon}->{'fsvars'};

			    
				    my $newclassesstr = $codons->{'pairs'}->{$bestcodon}->{'cluster_attrs'};
				    print "#CODON $bestcodon $codons->{'pairs'}->{$bestcodon}->{'gfreq'} max_annotated_len:$maxlen ";
				    if($debug){
					print "delta_len_max:$deltafracmax\n"; 
				    }
				    else{
					print "\n";
				    }
				    
				    #EDITTBL cluster_id, codon, genome_freq, annotated_freq, len, neworfs, overlaps, comments
				    print "#EDITTBL\tC$cluster_id\t$bestcodon\t$codons->{'pairs'}->{$bestcodon}->{'gfreq'}";
				    
				    if(scalar(keys%{$codons->{'pairs'}->{$bestcodon}->{'neworfs'}}) > 0){
					print "(N:",scalar(keys%{$codons->{'pairs'}->{$bestcodon}->{'neworfs'}}),")";
				    }
				    if(scalar(keys %{$codons->{'pairs'}->{$bestcodon}->{'fsvars'}}) > 0){
					print "(F:",scalar(keys%{$codons->{'pairs'}->{$bestcodon}->{'fsvars'}}),")";
				    }
				    print "\t$codons->{'pairs'}->{$bestcodon}->{'afreq'}\t$codons->{'pairs'}->{$bestcodon}->{'len'}\t";
				    print scalar(keys %{$codons->{'pairs'}->{$bestcodon}->{'olaps'}}),"\t";
				    if(scalar(keys %{$codons->{'pairs'}->{$bestcodon}->{'olaps'}})){
					print "#OVERLAP ";
					$altcodons->{$bestcodon}->{'isoverlap'}=1;
				    }
				    if($codonlength eq $maxlen){
					print "#MAXLENEDIT ";
					$altcodons->{$bestcodon}->{'maxlen'}=1;
				    }
				    if($codons->{'pairs'}->{$bestcodon}->{'gfreq'} eq (scalar(keys %$mappedorgs)+scalar(keys %$unmappedorgs))){
					print "#FCONSISTENT ";
					$altcodons->{$bestcodon}->{'fcon'}=1;
				    }
				    print "\n";
				    #Collapse all indels into runs
					foreach my $org (keys %{$codons->{'pairs'}->{$bestcodon}->{'fsvars'}}){
					
					# my @coords = sort {$a->[0] <=> $b->[0]} (@{$codons->{'pairs'}->{$bestcodon}->{'fsvars'}->{$org}});
# 					my @runs;
# 					my $indelstr1;
# 					my $indelstr2;
# 					my $last;
# 					my $start;
# 					my $end;
# 					for(my $i=0;$i<@coords;$i++){
# 					    #print join(',',@{$coords[$i]}),"\n";
# 					    if($i==0){
# 						$start=$coords[$i]->[0];
# 					    }
# 					    elsif(abs($last+1 - $coords[$i]->[0]) > 1){
# 						push @runs,[$start,$last,$indelstr1,$indelstr2];
# 						$indelstr1="";
# 						$indelstr2="";
# 						$start=$coords[$i]->[0];
# 					    }
# 					    $last=$coords[$i]->[0];
# 					    $indelstr1.=$coords[$i]->[1];
# 					    $indelstr2.=$coords[$i]->[2];
# 					}
# 					push @runs,[$start,$last,$indelstr1,$indelstr2];
# 					#Remove runs that are multiple of 3
# 					my @fsruns;
# 					foreach my $r (@runs){
# 					    die if(length($r->[2]) != length($r->[3]));
# 					    if(length($r->[2])%3!=0){
# 						push @fsruns,$r;
# 					    }
# 					}
					my @fsruns = @{$codons->{'pairs'}->{$bestcodon}->{'fsvars'}->{$org}};
					print "#FS $org ";
					foreach my $r (@fsruns){
					    print "[$r->[0]-$r->[1] $r->[2]:$r->[3]] $r->[4]";
					}
					print "\n";
				    }

				    foreach my $org (keys %{$codons->{'pairs'}->{$bestcodon}->{'neworfs'}}){
					my $fmin = $codons->{'pairs'}->{$bestcodon}->{'neworfs'}->{$org}->{'fmin'};
					my $fmax = $codons->{'pairs'}->{$bestcodon}->{'neworfs'}->{$org}->{'fmax'};
					my $orient = $codons->{'pairs'}->{$bestcodon}->{'neworfs'}->{$org}->{'orient'};
					my $olapgenes = &getFeaturesByInterval($atree,$org,$fmin,$fmax,$orient);
					if(scalar (keys %$olapgenes)>0){
					    foreach my $feat (keys %$olapgenes){
						print STDERR "Unexpected genes found ",join(",",@{$olapgenes->{$feat}}),"\n";
					    }
					    #die;
					}
					print "#NEWORF $org $fmin,$fmax,",($fmax-$fmin),",$orient\n";
					$new_orfs->{$org}++;
				    }
				}
			    }
			    else{
				#print STDERR "#WARNING Codon $bestcodon has 0 frequency\n";
			    }
			}
			close EFILE;
		    }
		}
	    }
	    #Print cluster
	    $classesallstr = join(';',sort {$a cmp $b} keys %{$cluster_attrs});
	    $classes_all->{$classesallstr}->{'ngenes'} += scalar(keys %$mappedgenes);
	    $classes_all->{$classesallstr}->{'nclusters'}++;
	    $classes_all->{$classesallstr}->{'new_orfs'}+= scalar(keys %$new_orfs);

	    &reportCluster($query,$mappedorgs,$mappedgenes,$unmappedorgs,$unmappedgenes,$feat_attrs,$cluster_attrs,$seq_attrs,$new_orfs);

	    $classesstr = join(';',sort {$a cmp $b} keys %{$cluster_attrs});
	    $clusters->{$cluster_id}->{'alts'} = $altcodons;
	    $clusters->{$cluster_id}->{'codons'} = $codons->{'pairs'};

	    $classes_sum->{$classesstr}->{'ngenes'} +=scalar(keys %$mappedgenes);
	    $classes_sum->{$classesstr}->{'nclusters'}++;
	    $classes_sum->{$classesstr}->{'new_orfs'}+= scalar(keys %$new_orfs);
	    $neworfcount+=scalar(keys %$new_orfs);

	    $validcluster++;

	    print "#VALID\tCLUSTER_$cluster_id\tNum_organisms=",scalar(keys %$mappedorgs)+1,
	    "\tNum_genes=",scalar(keys %$mappedgenes),"\n" if($debug);;

	    #For unmapped genes, save the best overlapping alignment
	    foreach my $feat_name (keys %$unmappedgenes){
		if(!exists $mapped->{$feat_name}){
		    if(exists $unmapped->{$feat_name} #first alignment encountered
		       || $unmappedgenes->{$feat_name}->{'cov'} > $unmapped->{$feat_name}->{'cov'}){ #better coverage
			$unmapped->{$feat_name}->{'cov'} = $unmappedgenes->{$feat_name}->{'cov'}/$features->{$feat_name}->[3]; #%coverage over gene length
			if($unmappedgenes->{$feat_name}->{'len'}){
			    $unmapped->{$feat_name}->{'pid'} = $unmappedgenes->{$feat_name}->{'pid'}/$unmappedgenes->{$feat_name}->{'len'}; #%id over aligned length
			}
			else{
			    $unmapped->{$feat_name}->{'pid'} = 0;
			}
			$unmapped->{$feat_name}->{'len'} = $unmappedgenes->{$feat_name}->{'len'};
			$unmapped->{$feat_name}->{'WGA_cluster'} = $cluster_id;
		    }
		    else{
			die if(exists $unmappedgenes->{$feat_name} && !exists $unmapped->{$feat_name});
		    }
		} 
	    }
	}	
	else{
	    die "Feature $query not mapped but marked so" if(exists $mapped->{$query});
	    #Cluster is a singleton, skip it or print for debugging
	    print "#SKIPPED\t$query\tWGA$cluster_id\tNum_organisms=",scalar(keys %$mappedorgs),
	    "\tNum_genes=",scalar(keys %$unmappedgenes),"\n" if($debug);
	    #This cluster was skipped because it does not pass coverage cutoffs
	    #Optionally print
	    if($printskipped){
		#print "#$query\tWGA$cluster_id\t$currorg\tcov:",$qcov/($fmax-$fmin),"\tid:1\tspan:$fmin-$fmax\tlen:",$fmax-$fmin,"\n";
		foreach my $organism (sort {$a cmp $b} keys %$unmappedorgs){
		    if(ref $unmappedorgs->{$organism} && exists $unmappedorgs->{$organism}->{'features'}){
			my($start,$end) = &getspan($unmappedgenes,keys %{$unmappedorgs->{$organism}->{'features'}});
			my @ogenes = sort {$features->{$a}->[1] <=> $features->{$b}->[1]} (keys %{$unmappedorgs->{$organism}->{'features'}});
			my @ocovs = map {sprintf("%.2f",$unmappedgenes->{$_}->{'cov'}/$features->{$_}->[3])} (@ogenes);
			my @oids  = map {sprintf("%.2f",$unmappedgenes->{$_}->{'pid'}/$unmappedgenes->{$_}->{'len'})} (@ogenes);
			
			print "#",join(',',@ogenes),
			"\tWGA$cluster_id",
			"\t$organism",
			"\tcov:",join(',',@ocovs),
			"\tid:",join(',',@oids),
			"\tspan:$start-$end len:",$end-$start,
			"\n" if($debug);
		    }
		}
	    }
	}

	#foreach my $organism (keys %$new_orfs){
	#my $orfidx=0;
	#foreach my $alt (@{$new_orfs->{$organism}}){
	#$neworfcount++;
	#}
	#}


	if($autocorrect){
	    #Auto-correct cluster
	    my @neworfs;
	    foreach my $organism (sort {$a cmp $b} keys %$mappedorgs){
		my @ogenes = sort {$features->{$a}->[1] <=> $features->{$b}->[1]} (keys %{$mappedorgs->{$organism}});
		my $classes;
		my $longestorf=0;		
		my $longestpairc=0;		
		foreach my $gene (@ogenes){
		    if(exists $feat_attrs->{$gene}){
			foreach my $c (sort {$a cmp $b} keys %{$feat_attrs->{$gene}}){
			    $classes->{$c}++;
			}
		    }
		    $longestorf = ($features->{$gene}->[3] > $longestorf) ? $features->{$gene}->[3] : $longestorf;
		    $longestpairc = ($feat_attrs->{$gene}->{'pairfreq'} > $longestpairc) ? $feat_attrs->{$gene}->{'pairfreq'} : $longestpairc;
		}
		##
		
		#my @attrs = sort {$a cmp $b} keys %$classes;
		if(exists  $seq_attrs->{$organism}){
		    #Report alternative start sites if they result in a longer ORF
		    my @alts;
		    my $orfidx=0;
		    foreach my $alt (@{$seq_attrs->{$organism}}){
			#Report alternative starts or possible frameshifts
			if($alt =~ /alt_start/){
			    print "#$alt\n" if($debug);;
			    #Only report if results in a longer ORF
			    my($astart,$aend,$aorient,$alen) = ($alt =~ /alt_start=(\d+)-(\d+),orient\:([^,]+),len\:(\d+)/);
			    my($apairfreq) = ($alt =~ /pairfreq:(\d+)/);
			    print "#alt $astart,$aend,$aorient,$alen\n" if($debug);;
			    print STDERR "BAD $alt" if(!$astart || !$aend || !$aorient || !$alen);
			    if(!$longer_altstarts || $alen>$longestorf){
				if(!$moreconsistent_altstarts || $apairfreq>$longestpairc){
				    push @alts,["ALTSTARTgene$organism$orfidx",$astart,$aend,$aend-$astart,$aorient];
				    $orfidx++;
				}
			    }
			    else{
				print "#Skipping $alt $alen<$longestorf\n" if($debug);;
			    }
			}
		    }
		    #Report frameshifts if they result in a longer ORF
		    #This should also include alt start, frameshift pairs if they result in a longer ORF
		    foreach my $alt (@{$seq_attrs->{$organism}}){
			#Report alternative starts or possible frameshifts
			if($alt =~ /alt_fs/){
			    print "#$alt\n" if($debug);;
			    #Only report if results in a longer ORF
			    my($astart,$aend,$aorient,$alen) = ($alt =~ /alt_fs=(\d+)-(\d+),orient\:([^,]+),len\:(\d+)/);
			    print "#alt $astart,$aend,$aorient,$alen\n" if($debug);;
			    die "$alt" if(!$astart || !$aend || !$aorient || !$alen);
			    if(!$longer_altstarts || $alen>$longestorf){
				push @alts,["ALTFSgene$organism$orfidx",$astart,$aend,$aend-$astart,$aorient];
				$orfidx++;
			    }
			    else{
				print "#Skipping $alt $alen<$longestorf\n" if($debug);;
			    }
			}
		    }
		    
		    
		    #Replace $ogenes
		    if(scalar(@alts)>0){
			my @sortedalts = sort {$b->[3] <=> $a->[3]} @alts;
			my $neworf = $sortedalts[0];
			print "#Num genes ",scalar(keys %$mappedgenes),"\n" if($debug);;
			if($autocorrect){
			    foreach my $gene (@ogenes){
				$deleted->{$gene}++;
				delete $mappedgenes->{$gene};
				delete $mappedorgs->{$organism}->{$gene};
				print "#Possible deleting $gene\n" if($debug);;
			    }
			}
			my $featlen = $neworf->[3];
			
			#if($neworfcov/$featlen >= $coverage_cutoff  && #%coverage over matching gene length
			#  $neworfpid/$alnlen >= $pid_cutoff){ #%id over aligned length onl
			$features->{$neworf->[0]} = [$organism,$neworf->[1],$neworf->[2],$neworf->[3],$neworf->[4]];
			$mappedgenes->{$neworf->[0]}->{'fmin'} = $neworf->[1];
			$mappedgenes->{$neworf->[0]}->{'fmax'} = $neworf->[2];
			$mappedgenes->{$neworf->[0]}->{'len'} = $neworf->[3];
			$mappedgenes->{$neworf->[0]}->{'relorient'} = $neworf->[4];
			$mappedorgs->{$organism}->{'features'}->{$neworf->[0]}++;
			$mappedorgs->{$organism}->{'qcov'} = '?';
			#Add new gene
			print "#Adding ",join(',',@$neworf),"\n" if($debug);;
			push @neworfs,$neworf->[0];
			$adjustedorfs++;
			print "#Num genes ",scalar(keys %$mappedgenes),"\n" if($debug);
		    }
		}
	    }
	    if(scalar(keys %$mappedgenes)>1 && scalar(keys %$mappedorgs)>1){
		my($feat_attrs,$cluster_attrs,$codons) = &annotateCluster($atree,$mappedgenes,$mappedorgs);
		my $classesstr = join(';',sort {$a cmp $b} keys %{$cluster_attrs});
		$newclasses_sum->{$classesstr}->{'ngenes'} +=scalar(keys %$mappedgenes);
		$newclasses_sum->{$classesstr}->{'nclusters'}++;
	    }
	    foreach my $neworf (@neworfs){
		delete $mappedgenes->{$neworf};
		delete $features->{$neworf};
	    }
	}
    }
}



print "#NUM CLUSTERS $validcluster\n";

#Mark the remaining features as singletons categorized as
#1)not found in any alignments !exists mapped && !exists unmapped
#2)aligned but below cutoffs   !exists mapped && exists unmapped
#If duplications file provided, mark accordingly
($nomatches,$dups) = &findSingletons($atree,$mapped,$unmapped,$subsumed,$datree);

#Calculate summary stats 
my $avgcov=0;
my $avgid=0;
my $mappedgenescount=0;

my $unmappedgenescount=0;
my $avgunmappedcov=0;
my $avgunmappedid=0;
my $unmappeddups=0;

my $nohit=0;
my $nohitdupcount=0;
my $neworfscount=0;

foreach my $feat_name (keys %$features){
    my $fmin = $features->{$feat_name}->[1];
    my $fmax = $features->{$feat_name}->[2];
    if(exists $mapped->{$feat_name}){
	die if(exists $unmapped->{$feat_name});
	die if(exists $nomatches->{$feat_name});
	if($mapped->{$feat_name}->{'cov'}>1){
	    print STDERR "Bad cov ",$mapped->{$feat_name}->{'cov'},"\n" if($verbose);
	    #$mapped->{$feat_name}->{'cov'}=1;
	}
	if($mapped->{$feat_name}->{'pid'}>1){
	    print STDERR "Bad id ",$mapped->{$feat_name}->{'pid'},"\n" if($verbose);
	    #$mapped->{$feat_name}->{'pid'}=1;
	}
	$avgcov+=$mapped->{$feat_name}->{'cov'};
	$avgid+=$mapped->{$feat_name}->{'pid'};
	$mappedgenescount++;
    }
    elsif(exists $unmapped->{$feat_name}){
	die if(exists $mapped->{$feat_name});
	die if(exists $nomatches->{$feat_name});	
	if($unmapped->{$feat_name}->{'cov'}>1){
	    print STDERR "Bad cov ",$unmapped->{$feat_name}->{'cov'},"\n" if($verbose);
	    #$unmapped->{$feat_name}->{'cov'}=1;
	}
	if($unmapped->{$feat_name}->{'pid'}>1){
	    print STDERR "Bad id ",$unmapped->{$feat_name}->{'pid'},"\n" if($verbose);
	    #$unmapped->{$feat_name}->{'pid'}=1;
	}
	
	$avgunmappedcov+=$unmapped->{$feat_name}->{'cov'};
	$avgunmappedid+=$unmapped->{$feat_name}->{'pid'};
	$unmappedgenescount++;
	if(exists $dups->{$feat_name}){
	    $unmappeddups++;
	}
    }
    elsif(exists $dups->{$feat_name}){
	$nohitdupcount++;
    }
    elsif(exists $nomatches->{$feat_name}){ 
	$nohit++;
    }
    else{
	#Genes should be categorized in mapped,unmapped,singletons
	die if(exists $mapped->{$feat_name});
	die if(exists $unmapped->{$feat_name});
	die if(exists $nomatches->{$feat_name});	
	#The rest are either deleted or newly called ORFs that are discarded
	die if(!exists $subsumed->{$feat_name} && !exists $neworfs->{$feat_name});
    }
}

print STDERR "#Mismatch between NOHIT=$nohit and nomatches lookup",scalar(keys %$nomatches),"\n" if($nohit != scalar(keys %$nomatches));
   
#Print summary stats
print "\n\n\n";
print "Class legend\n";
print "C{S,E}1 - consistent start,stop\n";
print "C{S,E}2 - inconsistent start,stop. More than one annotated in a group\n";
print "C{S,E}3 - unaligned start,stop\n";
print "C{S,E}4 - invalid start,stop according to translation table\n";
print "C{S,E}0 - start,stop at/off contig boundary\n";
print "CM1     - multiple gene fragments. possible interrupted genes\n";
print "CX      - invalid translation\n";
#print "CS/E2.1 - there is only one annotated gene for each genome, but not all genomes use the same start/stop\n";
#print "CS/E2.2 - the start/stop of some genes fall in a gapped region of the alignment\n";

print "Summary classes\n";
foreach my $cstr (sort {$classes_sum->{$b}->{'ngenes'} <=> $classes_sum->{$a}->{'ngenes'}} (keys %$classes_sum)){
    print "$cstr: num_genes:$classes_sum->{$cstr}->{'ngenes'} num_clusters:$classes_sum->{$cstr}->{'nclusters'}\n";
}

print "Complete classes\n";
foreach my $cstr (sort {$classes_all->{$b}->{'ngenes'} <=> $classes_all->{$a}->{'ngenes'}} (keys %$classes_all)){
    print "$cstr: num_genes:$classes_all->{$cstr}->{'ngenes'} num_clusters:$classes_all->{$cstr}->{'nclusters'}\n";
}
print "Number of clusters containing aligned features\n";
print "CLUSTERS:$validcluster\n";
print "Number aligned features mapped into clusters\n";
print "MAPPED:$mappedgenescount AVGCOV:",$avgcov/$mappedgenescount," AVGID:",$avgid/$mappedgenescount,"\n" if($mappedgenescount);
print "Number features with an overlapping alignment but are not mapped into clusters\n";
print "UNMAPPED:$unmappedgenescount AVGCOV:",$avgunmappedcov/$unmappedgenescount," AVGID:",$avgunmappedid/$mappedgenescount," NUMDUPS:$nohitdupcount\n" if($unmappedgenescount && $mappedgenescount);
if(exists $options{'duplications'}){
    print "Number of features with no mapping and marked as duplications\n";
    print "DUPS:$nohitdupcount\n";
}
print "Number of features with no overlapping alignment\n";
print "NOHIT:$nohit\n";
print "Number of missing annotations\n";
print "MISSORF:$neworfcount\n";

close $cfh;
close $ctfh;
close $ctfh2;

&printExtJS($clusters);

exit(0);

#############################
# Subroutines
#############################
#Primary method of obtaining mapped annotation from an alignment
#Build a cluster of aligned features/genes based on a single query
#gene, $query
#TODO: Confirm qcov,qpid,cov,pid are calculated correctly. Correct for overlapping alignments
sub buildCluster{
    my ($atree,$query) = @_;

    #Attributes of the query
    my $qseqname = $features->{$query}->[0];
    $feat2organism->{$query} = $qseqname;
    my $qcurrorg = '?';
    my $qcov = 0;
    my $qpid = 0;
    my $qalnfmin = undef;
    my $qalnfmax = undef;
    my $qfmin = $features->{$query}->[1];
    my $qfmax = $features->{$query}->[2];
    my $qfeatlen = $qfmax-$qfmin;
    my $qrelorient = 0;

    my $qorient = $features->{$query}->[4];
    
    print "#MAPFEATURE Mapping $query $qseqname:$qfmin-$qfmax len:",$qfmax-$qfmin,"\n" if($debug);;
    
    #AlignmentTree::map() 
    #returns [0=alignment_name,1=seqid,2=align_start,3=align_stop,4=align_cov,5=feature_name,6=seqid,7=feature_cov,8=feature_pid]

    my @isect = $atree->map($qseqname,$qfmin,$qfmax);
    
    #List of alignments that comprise the current cluster
    my $goodalignments = {};
    #List of seqs that overlap query in the cluster
    my $allseqs = {};

    #List of organism_ids in the current cluster 
    my $mappedorgs = {}; #passes cutoffs
    my $unmappedorgs = {}; #do not pass cutoffs
    
    #List of genes in the current cluster
    my $mappedgenes = {}; #passes cutoffs
    my $unmappedgenes = {}; #do not pass cutoffs
    
    #Contains list of annotations that are overlapping in an alignment
    my $alnfeats = {};
    my $alnorgs = {};

    my $valid=0;
    
    my $nisect;
    ($nisect,$allseqs,$goodalignments) = &getAlignedFeatures($atree,$qseqname,$query,$qfmin,$qfmax,'gene');
    

    if($verbose){
	print "#QUERY=$query len=$features->{$query}->[3]";
	print " coords=$qfmin-$qfmax len=$features->{$query}->[3] strand=$features->{$query}->[4]";
	print " Num_alignments=",scalar(keys %$goodalignments);
	print "\n";
    }
   
    #Report annotated frame relative to query
    my $seqalnpos;
    my @isect;
    if($qorient eq '+'){
	@isect = $atree->intersect($qseqname,$qfmin,$qfmin+3);
    }
    else{
	@isect = $atree->intersect($qseqname,$qfmax-3,$qfmax);
    }
    foreach my $r (@isect){
	my $feat_name = $r->[0];
	my $seqname = $r->[1];
	my $align_name = $r->[5];
	#print "#Setting $seqname query frame: qorient $qorient, alnframe: $r->[7] $r->[2]-$r->[3]\n";
	if(exists $goodalignments->{$feat_name} && $feat_name =~ /^WGA/){
	    if($qorient eq '+'){
		if($r->[7] eq '-'){
		    $seqalnpos->{$seqname}=$r->[3];
		}
		else{
		    $seqalnpos->{$seqname}=$r->[2];
		}
	    }
	    else{
		if($r->[7] eq '-'){
		    $seqalnpos->{$seqname}=$r->[3];
		}
		else{
		    $seqalnpos->{$seqname}=$r->[2];
		}
	    }
	}
    }

    foreach my $r ( sort { $features->{$b->[0]}->[3] <=> $features->{$a->[0]}->[3] } #sort on feature length
		    @$nisect){
	my $feat_name = $r->[0];
	my $seqname = $r->[1];
	my $align_name = $r->[5];
	#Check if we want to consider this alignment
	if(exists $goodalignments->{$align_name}){
	    my($alnobj,$bv,$width) = $atree->getAlignment($align_name);
	    $feat_name =~ s/gene\://;
	    if(!exists $features->{$feat_name}){
		print "#Bad feature found $feat_name. Not in input file. Skipping\n" if($debug);
		next;
	    }
	    #Capture some stats on the matching genes
	    #TODO the cov,pid stats assume non-overlapping alignments
	    if($query ne $feat_name){
		#Only report genes that have not been mapped
		if(!exists $mapped->{$feat_name} && !exists $deleted->{$feat_name} && exists $features->{$feat_name}){
		    print "#MAP:",join("\t",$cluster_id,@$r),"\n" if($debug);		
		    die "Mismatching orientation for $feat_name. Mapping showing $r->[12]. Input reporting $features->{$feat_name}->[4]" if($r->[12] ne $features->{$feat_name}->[4]);
		    die "fmax < fmin" if($r->[3]<$r->[2]);
		    die "Mismatched strand for $feat_name. Expecting $r->[12], got $features->{$feat_name}->[4]" if($r->[12] ne $features->{$feat_name}->[4]);
		    #Sum the coverage for each gene versus the query
		    if(exists $alnfeats->{$feat_name}->{'fmin'}){
			$alnfeats->{$feat_name}->{'fmin'}=($r->[2]<$alnfeats->{$feat_name}->{'fmin'}) ? $r->[2]: $alnfeats->{$feat_name}->{'fmin'};
		    }
		    else{
			$alnfeats->{$feat_name}->{'fmin'}=$r->[2];
		    }
		    if(exists $alnfeats->{$feat_name}->{'fmax'}){
			$alnfeats->{$feat_name}->{'fmax'}=($r->[3]>$alnfeats->{$feat_name}->{'fmax'}) ? $r->[3]: $alnfeats->{$feat_name}->{'fmax'};
		    }
		    else{
			$alnfeats->{$feat_name}->{'fmax'}=$r->[3];
		    }
		    $alnfeats->{$feat_name}->{'cov'}+=$r->[7];
		    $alnfeats->{$feat_name}->{'pid'}+=$r->[8];
		    $alnfeats->{$feat_name}->{'len'}+=($r->[3]-$r->[2]);
		    die "Bad pid $alnfeats->{$feat_name}->{'pid'} > $alnfeats->{$feat_name}->{'len'} from pid:$r->[8] len:($r->[3]-$r->[2]) ".($r->[3]-$r->[2]) if($alnfeats->{$feat_name}->{'pid'} > $alnfeats->{$feat_name}->{'len'});
		    $alnfeats->{$feat_name}->{'relorient'} = $r->[11]; 
		    $feat2organism->{$feat_name} = $r->[1];
		    #num aligned residues $r->[8] indicates matches on query seq
		    #   |NNNNNN---NNNNNNNNN| query 15 residues
		    #      |NNNNNNNNN---NNNNNN| hit 15 residues - 9 matching , qcov=9/15, cov=9/15
		    $alnfeats->{$feat_name}->{'qcov'}+=$r->[4];
		    print "#$feat_name $align_name $r->[2]-$r->[3] len:",$r->[3]-$r->[2],
		    ",$alnfeats->{$feat_name}->{'len'} cov:$r->[7],$alnfeats->{$feat_name}->{'cov'} id:$r->[8],$alnfeats->{$feat_name}->{'pid'} alnorient:$r->[10] featorient:$r->[11]\n" if($debug);
		}
		else{
		    #This feature has already been mapped
		    print "#Alternative mapping for $feat_name cov:$r->[4] pid:$r->[8] len:",$r->[3]-$r->[2]," matchingorient:$r->[10],$r->[11]\n" if($debug);;
		}
	    }
	    else{
		die if($feat_name ne $query);
		die if($r->[10] ne $r->[11]);
		#Capture some stats on the query
		if(defined $qalnfmin){
		    $qalnfmin = ($r->[2] < $qalnfmin) ? $r->[2] : $qalnfmin;
		}
		else{
		    $qalnfmin = $r->[2];
		}
		if(defined $qalnfmax){
		    $qalnfmax = ($r->[3] > $qalnfmax) ? $r->[3] : $qalnfmax;
		}
		else{
		    $qalnfmax = $r->[3];
		}
		$qcov += $r->[7];
		$qcurrorg = $r->[6];
		$qpid += $r->[8];
		$qrelorient = $r->[10];
	    }
	}
    }

    $mappedgenes->{$query}->{'fmin'} = $qalnfmin;
    $mappedgenes->{$query}->{'fmax'} = $qalnfmax;
    $mappedgenes->{$query}->{'cov'} = $qalnfmax-$qalnfmin;#$qcov;
    $mappedgenes->{$query}->{'pid'} = $qpid;#TODO, not pid, rather %aln: allows mismatches but no gaps
    $mappedgenes->{$query}->{'len'} = $qfeatlen;
    $mappedgenes->{$query}->{'relorient'} = $qrelorient;
    $mappedgenes->{$query}->{'alignments'} = [keys %$goodalignments];
    $mappedorgs->{$qcurrorg}->{'features'}->{$query}++;
    $mappedorgs->{$qcurrorg}->{'qcov'} = $qalnfmax-$qalnfmin;


    #Set query coverage
    foreach my $feat_name (keys %$alnfeats){
	my $seqname = $features->{$feat_name}->[0];
	my $fmin = $features->{$feat_name}->[1];
	my $fmax = $features->{$feat_name}->[2];
	my $orient = $features->{$feat_name}->[4];
	my @isect;
	my $qmatchstart;
	if($orient eq '+'){
	    @isect = $atree->intersect($seqname,$fmin,$fmin+3);
	}
	else{
	    @isect = $atree->intersect($seqname,$fmax-3,$fmax);
	}
	my $alignedstart=0;
	foreach my $r (@isect){
	    my $feat_name = $r->[0];
	    if(exists $goodalignments->{$feat_name} && $feat_name =~ /^WGA/){
		if($r->[1] eq $qseqname){
		    my $align_name = $r->[5];
		    if($r->[7] eq '-'){
			$qmatchstart=$r->[3];
		    }
		    else{
			$qmatchstart=$r->[2];
		    }
		}
		elsif($r->[1] eq $seqname){
		    if(($orient eq '+' && $fmin == $r->[2] && $fmin+3 == $r->[3]) || 
		       ($orient eq '-' && $fmax-3 == $r->[2] && $fmax == $r->[3])){
			$alignedstart=1;
			print "#Aligned start $feat_name $seqname $orient $fmin-$fmax $r->[2]-$r->[3]\n" if($debug);
		    }
		    else{
			print "#Unaligned start for $feat_name $seqname $fmin == $r->[2] && $fmax == $r->[3]\n" if($debug);
		    }
		}
	    }
	}
	    
	my $featlen = $fmax-$fmin;
	if($alnfeats->{$feat_name}->{'cov'}/$featlen >= $coverage_cutoff  && #%coverage over matching gene length
	   $alnfeats->{$feat_name}->{'pid'}/$alnfeats->{$feat_name}->{'len'} >= $pid_cutoff){ #%id over aligned length only
	    print "Summing query coverage feat_name $feat_name $feat2organism->{$feat_name} = $alnfeats->{$feat_name}->{'qcov'}. Current total $alnorgs->{$feat2organism->{$feat_name}}->{'qcov'}\n" if($debug);
	    $alnorgs->{$feat2organism->{$feat_name}}->{'qcov'} += $alnfeats->{$feat_name}->{'qcov'};
	}

	if(exists $seqalnpos->{$seqname} && $alignedstart){
	    if($features->{$feat_name}->[4] eq '+' ){
		#print "#Feat $feat_name $features->{$feat_name}->[4] $seqalnpos->{$seqname}-$fmin ",$seqalnpos->{$seqname}-$fmin," ",($seqalnpos->{$seqname}-$fmax)%3,"\n";
		$alnfeats->{$feat_name}->{'relqrysdist'}=($seqalnpos->{$seqname}-$fmin);
		$alnfeats->{$feat_name}->{'frame'}=($seqalnpos->{$seqname}-$fmin)%3;
		$alnfeats->{$feat_name}->{'frameinqry'}=($qfmin-$qmatchstart)%3;
		#print "#qframe $feat_name $alnfeats->{$feat_name}->{'frameinqry'}\n";
	    }
	    else{
		#print "#Feat $feat_name $features->{$feat_name}->[4] $seqalnpos->{$seqname}-$fmax ",$seqalnpos->{$seqname}-$fmax," ",($seqalnpos->{$seqname}-$fmax)%3,"\n";
		$alnfeats->{$feat_name}->{'relqrysdist'}=($seqalnpos->{$seqname}-$fmax);
		$alnfeats->{$feat_name}->{'frame'}=($seqalnpos->{$seqname}-$fmax)%3;
		$alnfeats->{$feat_name}->{'frameinqry'}=($qfmax-$qmatchstart)%3;
		#print "#qframe $feat_name $alnfeats->{$feat_name}->{'frameinqry'}\n";
	    }
	}
    }
    
    foreach my $feat_name (keys %$alnfeats){
	#Check gene is part of input feature list [optional]
	die "Bad gene $feat_name" if(! exists $features->{$feat_name});
	die "Query gene should not map to itself" if($feat_name eq $query);
	#die "Can't find $feat_name in organism lookup" if(!exists $feat2organism->{$feat_name});
	#die "Bad organism $feat2organism->{$feat_name" if(!exists $alnorgs->{$feat2organism->{$feat_name}});
	#Check gene has not already been mapped
	die if(exists $mapped->{$feat_name});
	#coverage cutoff and percent identity cutoff
	my $fmin = $features->{$feat_name}->[1];
	my $fmax = $features->{$feat_name}->[2];
	my $featlen = $fmax-$fmin;
	die if($featlen<1);
	if($verbose){
	    print STDERR "Bad query coverage $qcov > $qfeatlen for $feat_name $fmin-$fmax\n" if($qcov > $qfeatlen);
	    print STDERR "Bad match coverage $alnfeats->{$feat_name}->{'cov'} > $featlen==$fmax-$fmin for $feat_name\n" if($alnfeats->{$feat_name}->{'cov'} > $featlen);
	    print STDERR "Bad match pid $alnfeats->{$feat_name}->{'pid'} > $alnfeats->{$feat_name}->{'len'} for $feat_name\n" if($alnfeats->{$feat_name}->{'pid'} > $alnfeats->{$feat_name}->{'len'});
	    print STDERR "#WARNING Bad len $alnfeats->{$feat_name}->{'len'} > ($features->{$feat_name}->[2]-$features->{$feat_name}->[1]) ".($features->{$feat_name}->[2]-$features->{$feat_name}->[1])." for $feat_name\n"
		if($alnfeats->{$feat_name}->{'len'} > ($features->{$feat_name}->[2]-$features->{$feat_name}->[1]));
	}
	#
	#Coverage and percent_id cutoffs are checked here in the following order 
	#Check that coverage over shorter of query and hit
	#query_coverage > coverage_cutoff || hit_coverage > coverage_cutoff && hit_pid > pid_cutoff
	#query_coverage > coverage_cutoff && hit_coverage > coverage_cutoff && hit_pid > pid_cutoff
	#
	print "Cutoff check $feat_name $feat2organism->{$feat_name} $alnorgs->{$feat2organism->{$feat_name}}->{'qcov'},$qfeatlen qcov=",($alnorgs->{$feat2organism->{$feat_name}}->{'qcov'}/$qfeatlen)," >= $query_coverage_cutoff ",
	$alnfeats->{$feat_name}->{'cov'}/$featlen ," >=  $coverage_cutoff ", 
	$alnfeats->{$feat_name}->{'pid'}/$alnfeats->{$feat_name}->{'len'},">= $pid_cutoff\n" if($debug);

	if(($query_coverage_cutoff==0 || ($alnorgs->{$feat2organism->{$feat_name}}->{'qcov'}/$qfeatlen >= $query_coverage_cutoff)) && #%coverage over query(longer feature in the comparison)
	   $alnfeats->{$feat_name}->{'cov'}/$featlen >= $coverage_cutoff  && #%coverage over matching gene length (shorter feature in the comparison)
	   $alnfeats->{$feat_name}->{'pid'}/$alnfeats->{$feat_name}->{'len'} >= $pid_cutoff){ #%id over aligned length only
	    
	    print "PASSED $feat_name\n" if($debug);
	    
	    #Check matching len is <= length of gene
	    $mappedorgs->{$feat2organism->{$feat_name}}->{'features'}->{$feat_name}++;
	    #print "WARNING query coverage > query length: $alnorgs->{$feat2organism->{$feat_name}}->{'qcov'} > $qfeatlen\n" 
	    #if($alnorgs->{$feat2organism->{$feat_name}}->{'qcov'} > $qfeatlen);
	    $mappedorgs->{$feat2organism->{$feat_name}}->{'qcov'} = $alnorgs->{$feat2organism->{$feat_name}}->{'qcov'};

	    $mappedgenes->{$feat_name}->{'fmin'} = $alnfeats->{$feat_name}->{'fmin'};
	    $mappedgenes->{$feat_name}->{'fmax'} = $alnfeats->{$feat_name}->{'fmax'};
	    $mappedgenes->{$feat_name}->{'cov'} = $alnfeats->{$feat_name}->{'cov'};
	    $mappedgenes->{$feat_name}->{'pid'} = $alnfeats->{$feat_name}->{'pid'};
	    $mappedgenes->{$feat_name}->{'len'} = $alnfeats->{$feat_name}->{'len'};
	    $mappedgenes->{$feat_name}->{'relorient'} = $alnfeats->{$feat_name}->{'relorient'};
	    #print "FRAME $alnfeats->{$feat_name}->{'frame'}\n";
	    #print "RELQRYSTARTDIST $alnfeats->{$feat_name}->{'relqrysdist'}\n";
	    $mappedgenes->{$feat_name}->{'frame'} = $alnfeats->{$feat_name}->{'frame'};
	    $mappedgenes->{$feat_name}->{'frameinqry'} = $alnfeats->{$feat_name}->{'frameinqry'};
	    $mappedgenes->{$feat_name}->{'relqrysdist'} = $alnfeats->{$feat_name}->{'relqrysdist'};
	}
	else{
	    print "BELOW $feat_name\n" if($debug);
	    #Does not pass cutoffs	    
	    $unmappedgenes->{$feat_name}->{'cov'} = $alnfeats->{$feat_name}->{'cov'};
	    $unmappedgenes->{$feat_name}->{'fmin'} = $alnfeats->{$feat_name}->{'fmin'};
	    $unmappedgenes->{$feat_name}->{'fmax'} = $alnfeats->{$feat_name}->{'fmax'};
	    $unmappedgenes->{$feat_name}->{'pid'} = $alnfeats->{$feat_name}->{'pid'};
	    $unmappedgenes->{$feat_name}->{'len'} = $alnfeats->{$feat_name}->{'len'};
	    $unmappedgenes->{$feat_name}->{'relorient'} = $alnfeats->{$feat_name}->{'relorient'};
	    $unmappedgenes->{$feat_name}->{'frame'} = $alnfeats->{$feat_name}->{'frame'};
	    $unmappedgenes->{$feat_name}->{'frameinqry'} = $alnfeats->{$feat_name}->{'frameinqry'};
	    $unmappedgenes->{$feat_name}->{'relqrysdist'} = $alnfeats->{$feat_name}->{'relqrydist'};
	}
    }
    foreach my $seq (keys %$allseqs){
	if(!exists $mappedorgs->{$seq}){
	    #Does not pass cutoffs
	    $unmappedorgs->{$seq} = {};
	}
    }
    foreach my $feat_name (keys %{$unmappedgenes}){
	if(!exists $mappedorgs->{$feat2organism->{$feat_name}}){
	    die "ORG found in mapped list $feat2organism->{$feat_name} $feat_name query:$query queryorg:$qcurrorg" if(exists $mappedorgs->{$feat2organism->{$feat_name}});
	    $unmappedorgs->{$feat2organism->{$feat_name}}->{'features'}->{$feat_name}++;
	    $unmappedorgs->{$feat2organism->{$feat_name}}->{'qcov'} += $alnfeats->{$feat_name}->{'qcov'};
	}
    }
    return($mappedorgs,$mappedgenes,$unmappedorgs,$unmappedgenes);
}
    


#Classify consistency of annotations within a cluster
#Returns start,stop codon positions of annotated genes only
#
#Clusters are assigned one or more classes based on consistent gene structures
#Class CS1: All start codons in the cluster are aligned
#Class CS2: There are multiple, inconsistent start codons in the cluster
#Class CS3: One or more of the start codons are not aligned in the cluster
#Class CS4: Invalid annotated start codon
#Class CE1-3. Same as CS1-3 but for stop codons
#Class CM1 : Multiple spanned features in the cluster
sub annotateCluster{
    my($atree,$genes,$orgs) = @_;

    my $cluster_attrs = {};
    my $feat_attrs = {};

    my $starts = {};
    my $stops = {};
    my $codonpairs = {};

    my $alignedstartcount=0;
    my $alignedstopcount=0;

    my $seqstarts = {};
    my $seqstops = {};
    my $featstarts = {};
    my $featstops = {};

    foreach my $org (keys %$orgs){
	if(scalar(keys %{$orgs->{$org}->{'features'}})>1){
	    print "#Class CM1. Multiple genes spanning query. Count ",scalar(keys %{$orgs->{$org}->{'features'}}),"\n" if($debug);;
	    $cluster_attrs->{'CM1'} = [$org,scalar(keys %{$orgs->{$org}->{'features'}})];
	}
    }
    print "#Annotating cluster\n" if($debug);;
    foreach my $feat_name (keys %$genes){
	#Save relative position of start and stop codons in the
        #alignment $align_name
	die if(!exists $features->{$feat_name});
	my ($seqname,$fmin,$fmax,$len,$orient) = @{$features->{$feat_name}};
	my $relorient = $genes->{$feat_name}; #relative orientation of the annotation on the aligned seq
	#$relorient == 1 Annotation and alignment are on the same strand
	#$relorient == 0 Annotation and alignment are on opposite strands
	my($startcodon,$stopcodon,$partial_start,$partial_stop,$bad_start,$bad_stop) = &findCodons($atree,
												   $seqname,
												   $fmin,
												   $fmax,
												   $orient,$feat_name);
	if($verbose && !$bad_stop && !$bad_start){
	    print "BAD ORF $seqname,$fmin,$fmax\n" if(&isORF($db,$seqname,$fmin,$fmax,$orient)<=0);
	}


	if(ref $startcodon){
	    my($mcol,$align_name) = (@$startcodon);
	    my $token = $mcol.$CODON_DELIM.$align_name;
	    if($debug){
		if($orient eq '+'){
		    my @res= AlignmentTree::coordstocolumn($atree->{_alignments}->{$align_name}->[0],$seqname,$fmin,$fmin+3);
		    die "$res[0] ne $mcol $seqname,$fmin,$fmin+3" if($res[0] ne $mcol);
		}
		else{
		    my @res= AlignmentTree::coordstocolumn($atree->{_alignments}->{$align_name}->[0],$seqname,$fmax-3,$fmax);
		    die "$res[0] ne $mcol $seqname,$fmin,$fmin+3" if($res[0] ne $mcol);
		}
	    }
	    $starts->{$token}++;
	    print "#Start codon $feat_name $startcodon->[0] $startcodon->[1] $startcodon->[2] $startcodon->[3] $orient\n" if($debug);;
	    $alignedstartcount++;
	    $features->{$feat_name}->[7] = $startcodon->[0];
	    $features->{$feat_name}->[8] = $startcodon->[1];
	    $seqstarts->{$seqname}->{$token}++;
	    $featstarts->{$feat_name} = $token;
	    if($partial_start){
		$feat_attrs->{$feat_name}->{'CS0'}++; #start codon in PMARK spacer adjacent to contig boundary
		$cluster_attrs->{'CS0'}++; 
	    }
	    if($debug){
		$feat_attrs->{$feat_name}->{'startcol:'.$mcol}++;
	    }
	    if($bad_start){#$startcodon == -1){
		$feat_attrs->{$feat_name}->{'CS4'}++; #invalid start
		$cluster_attrs->{'CS4'}++; 
	    }

	}
	else{
	    $feat_attrs->{$feat_name}->{'CS3'}++;
	    $cluster_attrs->{'CS3'}++; 
	}
	if(ref $stopcodon){
	    my($mcol,$align_name) = (@$stopcodon);
	    my $token = $mcol.$CODON_DELIM.$align_name;
	    if($debug){
		if($orient eq '+'){
		    my @res= AlignmentTree::coordstocolumn($atree->{_alignments}->{$align_name}->[0],$seqname,$fmax-3,$fmax);
		    die "$res[0] ne $mcol" if($res[0] ne $mcol);
		}
		else{
		    my @res= AlignmentTree::coordstocolumn($atree->{_alignments}->{$align_name}->[0],$seqname,$fmin,$fmin+3);
		    die "$res[0] ne $mcol" if($res[0] ne $mcol);
		}
	    }
	    $stops->{$token}++;
	    print "#Stop  codon $feat_name $stopcodon->[0] $stopcodon->[1] $stopcodon->[2] $stopcodon->[3] $orient\n" if($debug);;
	    $alignedstopcount++;
	    $features->{$feat_name}->[9] = $stopcodon->[0];
	    $features->{$feat_name}->[10] = $stopcodon->[1];
	    $seqstops->{$seqname}->{$token}++;
	    $featstops->{$feat_name} = $token;
	    if($partial_stop){
		$feat_attrs->{$feat_name}->{'CE0'}++; #stop codon in PMARK spacer adjacent to contig boundary
		$cluster_attrs->{'CE0'}++; 
	    }
	    if($debug){
		$feat_attrs->{$feat_name}->{'stopcol:'.$mcol}++;
	    }
	    if($bad_stop){#$stopcodon == -1){
		$feat_attrs->{$feat_name}->{'CE4'}++; #invalid stop
		$cluster_attrs->{'CE4'}++; 
	    }
	}
	else{
	    $feat_attrs->{$feat_name}->{'CE3'}++;
	    $cluster_attrs->{'CE3'}++; 
	}
	if(scalar(keys%{$orgs->{$seqname}->{'features'}}==1) && exists $featstarts->{$feat_name} && $featstops->{$feat_name}){
	    #$codonpairs->{$featstarts->{$feat_name}.':'.$featstops->{$feat_name}}->{'gfreq'}++;
	    #$codonpairs->{$featstarts->{$feat_name}.':'.$featstops->{$feat_name}}->{'afreq'}++;
	    #$codonpairs->{$featstarts->{$feat_name}.':'.$featstops->{$feat_name}}->{'length'}+=$len;

	    $codonpairs->{$featstarts->{$feat_name}.':'.$featstops->{$feat_name}}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,1]; #[fmin,fmax,orient,is_annotated,fs_type]
	    
	}

    }

    if(scalar(keys %$starts)==1){
	#There is only one annotated start
	my @start = keys %$starts; 
	if($starts->{$start[0]}==scalar(keys %$genes)){
	    #and every gene has this annotated start
	    print "#Class CS1. Consistent starts\n" if($debug);;
	    $cluster_attrs->{'CS1'}++;
	}
	else{
	    #some genes are missing this start codon but there are no others
	    print "#Class CS3. Unaligned starts ",$starts->{$start[0]}, "==",scalar(keys %$genes),"\n" if($debug);;
	    $cluster_attrs->{'CS2'}++; 
	    $cluster_attrs->{'CS3'}++; 
	}
    }
    else{
	if($alignedstartcount == scalar(keys %$genes)){
	    #there is one annotated start codon for each genome, but not all genomes use the same start
	    print "#Class CS2. Inconsistent starts\n" if($debug);;
	    $cluster_attrs->{'CS2'}++;
	}
	else{
	    #there are multiple annotated start codons for genome
	    print "#Class CS3. Unaligned starts ",$alignedstartcount," == ",scalar(keys %$genes),"\n" if($debug);;
	    $cluster_attrs->{'CS2'}++; 
	}
    }
    if(scalar(keys %$stops)==1){
	#There is only one annotated stop
	my @stop = keys %$stops;
	if($stops->{$stop[0]}==scalar(keys %$genes)){
	    #and every gene is annotated with this stop
	    print "#Class CE1. Consistent stops\n" if($debug);;
	    $cluster_attrs->{'CE1'}++;
	}
	else{
            #some genes are missing this stop codon but there are no others
	    print "#Class CE3. Unaligned stops\n" if($debug);;
	    $cluster_attrs->{'CE2'}++; 
	    $cluster_attrs->{'CE3'}++; 
	}
    }
    else{
	if($alignedstopcount == scalar(keys %$genes)){
	    #there is one annotated stop codon for each genome, but not all genomes use the same stop
	    print "#Class CE2. Inconsistent stops\n" if($debug);;
	    $cluster_attrs->{'CE2'}++;
	}
	else{
	    #there are multiple annotated stop codons for genome
	    print "#Class CE3. Unaligned stops\n" if($debug);;
	    $cluster_attrs->{'CE2'}++; 
	}
    }
    
    #Save frequency of annotated starts, stops
    foreach my $feat_name (keys %$genes){
	#$feat_attrs->{$feat_name}->{'pairfreq='.$codonpairs->{"$featstarts->{$feat_name}"."$featstops->{$feat_name}"}}++;
	if(exists $featstarts->{$feat_name}){
	    $feat_attrs->{$feat_name}->{'startfreq='.$starts->{$featstarts->{$feat_name}}}++;
	    $feat_attrs->{$feat_name}->{'startcodon='.$featstarts->{$feat_name}}++;
	}
	if(exists $featstops->{$feat_name}){
	    $feat_attrs->{$feat_name}->{'stopfreq='.$stops->{$featstops->{$feat_name}}}++;
	    $feat_attrs->{$feat_name}->{'stopcodon='.$featstops->{$feat_name}}++;
	}
    }
    return ($feat_attrs,$cluster_attrs,{'starts'=>$seqstarts,'stops'=>$seqstops,'pairs'=>$codonpairs,'featstops'=>$featstops,'featstarts'=>$featstarts});
}

###################################
#Classify singletons and unannotated regions
#
#Singletons consist of all annotated ORFs that do not map into an existing cluster above cutoffs
#Singletons are classified into the following classes
#Class SLTN1: there are no alignments that overlap the singleton. apparently true singleton
#Class SLTN2: there are overlapping alignments, annotated ORF start can be modified to pass cutoffs into an existing cluster
#Class SLTN3: there are overlapping alignments, annotated ORF stop can be modified to pass cutoffs into an existing cluster
#Class SLTN4: there are overlapping alignments and unannotated ORFs can be mapped above cutoffs
#Class SLTN5: there are overlapping alignments, but no overlapping ORFs above cutoffs
sub annotateSingletons{
    my($atree,$seqname,$feat_name,$fmin,$fmax) = @_;
    my @classes;
    my @isect = $atree->intersect($seqname,$fmin,$fmax,$aligntoken);
    my $goodalignments = {};
    foreach my $r (@isect){
	my $feat_name = $r->[0];
	my $seqname = $r->[1];
	my $align_name = $r->[5];
	#Only consider WGA alignments (alignment name in $align_name) that span query (gene name in $feat_name)
	if($feat_name eq 'gene:'.$feat_name){
	    $goodalignments->{$align_name}++;
	}
    }
    if(scalar (@isect)==0){
	push @classes,"classSLTN1";
    }
    else{
	push @classes,"classSLTN5 Num_alns:".scalar(keys %$goodalignments);
    }
    return \@classes;
}

#Check if fmin-fmax,orient on seqname is a valid ORF
sub isORF(){
    my($db,$seqname,$fmin,$fmax,$orient,$fs) = @_;
    #hack to avoid some bioperl warnings that i cannot turn off
    Bio::Root::Root::verbose(0);
    *TEMP = \*STDERR;
    open(FOO, ">/dev/null");
    my $seqlen = ($fmax-$fmin);
    if((! defined $fs && $seqlen%3!=0) || $seqlen > $MAXORFLEN || $seqlen < $MINORFLEN){
	print "#Bad ORF length $seqname $fmin-$fmax ",$seqlen," ",$seqlen%3,"\n" if($verbose);
	return 0;
    }
    #my $seqobj = $db->get_Seq_by_id($seqname);
    my $seqobj = $db->{$seqname};
    die "Bad coordinates $fmin-$fmax @_" if($fmin >= $fmax);
    my $codon_table = Bio::Tools::CodonTable->new(-id=>11);
    if($seqobj){
	if($orient eq '+'){
	    #die "Bad coordinates $fmax extends past end of sequence" if($fmax >= $seqobj->length());

	    my $newobj;
	    my $fsadj=0;
	    my $pmark=0;
	    my $adj=0;
	    if($fs){
		my $newobjs = $seqobj->trunc($fmin+1,$fmax);
		my $gseq = $newobjs->seq();
		#print "Seq size ",length($gseq)," ",$newobjs->length(),"\n";
		#Check for PMARK spacer
		#my $encoding = 'C'x$newobjs->length();
		print "GSEQPRE:$gseq\n" if($debug);
		foreach my $f (sort {$b->[0] <=> $a->[0]} @$fs){
		    #print "SIZE ",scalar(@$fs),"\n";
		    foreach my $start (sort {$b <=> $a} @{$f->[5]}){
			if($start>=$fmin && $start<=$fmax){
			    #print "SAM $seqname,$fsadj $fmin-$fmax $start $f->[0] $f->[1] $f->[2] $f->[3] $f->[4] $adj\n";
			    #$fsadj+=$f->[4];
			    #die if(($start-$fmin) >= length($encoding));

			    die if(($start-$fmin) < 0);
			    
			    if($f->[4] == 1){
				#substr($encoding,$start-$fmin+$adj,1,'F');
				substr($gseq,$start-$fmin,1) = '';
			    }
			    elsif($f->[4] == -1){
				#substr($encoding,$start-$fmin+1+$adj,0) = 'B';
				substr($gseq,$start-$fmin+1,0) = 'N';
				$adj++;
			    }
			    elsif($f->[4] == 0){
				#substr($encoding,$start-$fmin,1,'G');
			    }
			}
		    }
		}

		if($newobjs->seq() =~ /$PMARK_SPACER/){
		    my $sloc = $-[0];
		    print "FOUND PMARK+ $sloc ",substr($newobjs->seq(),$sloc,36),"\n" if($verbose);
		#    substr($encoding,$sloc,36,'G'x36);
		    $pmark=1;
		}
		$newobj = new Bio::Seq(-seq=>$gseq);
		#$newobj = $newobj->revcom();
		print "GSEQPOST:$gseq\n" if($debug);
		#print "GSEQobj:",$newobj->seq(),"\n";
		#print "#Encoding $encoding\n" if($debug);
		#return 0 if(($seqlen+$fsadj)%3!=0);
		# if(0){
# 		    $newobjs->verbose(0);
# 		    Bio::Root::Root::verbose(0);
# 		    eval{
# 			$newobj = new Bio::Seq::EncodedSeq(-seq=>$newobjs->seq(),
# 							   -encoding=>$encoding,
# 							   -verbose=>0,
# 							   );
			
			
# 		    }
# 		    or do{
# 			print "ERROR: ",$@,"\n" if($verbose);
# 			print "$seqname $fmin,$fmax,$orient PMARK=$pmark ",$newobjs->seq(),"\n",$encoding,"\n" if($verbose);
# 			return 0;
			
# 		    };
# 		}
	    }
	    else{
		$newobj = $seqobj->trunc($fmin+1,$fmax);
	    }
	    
	    die if($verbose && $newobj->length() > $MAXORFLEN);
	    die if($verbose && $newobj->length() < $MINORFLEN);
	    
	    if(1){#if($codon_table->is_start_codon($newobj->subseq(1,3)) && 
		      #($codon_table->is_ter_codon($newobj->subseq($newobj->length()-3+1,$newobj->length())))){
		      #*STDERR = *FOO;
		my $protein_seq_obj;
		eval{
		    if(0 && $fs){
			#print "Using FS\n";
			$protein_seq_obj = $newobj->cds()->translate(
								     -codontable_id =>11,
								     #-orf=>1,
								     -complete => 1,
								     -throw => 1,
								     -verbose => 0
								     );
		    }
		    else{
			$protein_seq_obj = $newobj->translate(-codontable_id =>11,
							      -complete => 1,
							      -throw => 1
							      );
		    }
		}
		or do {
		    print "ERROR translate: ",$@,"\n" if($verbose);
		    print "ERROR translate $seqname $fmin,$fmax,$orient ",$newobj->seq(),"\n" if($verbose);
		    return 0;
		};
		#*STDERR = *TEMP;
		return 0 if(!$protein_seq_obj);
		if($protein_seq_obj->length()>0){# +1 == ($seqlen+$fsadj)/3){
		    die if($protein_seq_obj->seq() =~ /\*/);
		    return 1+$pmark;
		}
		else{
		    print "#Unexpected sequence length ",$protein_seq_obj->length()," expecting ",($seqlen+$fsadj)/3," from ORF $seqname $fmin-$fmax $orient ",$protein_seq_obj->seq(),"\n" if($verbose);
		}
	    }
	    else{
		#print "Possible alternative ORF on $seqname $fmin-$fmax,$orient has invalid start:",$newobj->subseq(1,3)," ",$codon_table->is_start_codon($newobj->subseq(1,3))," or stop:",$newobj->subseq($newobj->length()-3+1,$newobj->length())," ",$codon_table->is_ter_codon($newobj->subseq($newobj->length()-3+1,$newobj->length())),"\n" if($verbose);
		return -1;
	    }
	}
	else{
	    die if($orient ne '-');
	    my $newobj;
	    my $fsadj=0;
	    my $pmark=0;
	    my $adj=0;
	    my $encoding;
	    if($fs){
		my $newobjs = $seqobj->trunc($fmin+1,$fmax);
		my $gseq = $newobjs->seq();
		#$encoding = 'C'x$newobjs->length();
		print "GSEQPRE:$gseq\n" if($debug);
		#print "FS ",join(',',@$fs),"\n";
		#print $newobjs->seq(),"\n";
 		foreach my $f (sort {$b->[0] <=> $a->[0]} @$fs){
		    #print "SIZE ",scalar(@$fs),"\n";
		    foreach my $start (sort {$b <=> $a} @{$f->[5]}){
			if($start>=$fmin && $start<=$fmax){
			    #print "SAM $seqname,$fsadj $fmin-$fmax $start $f->[4] $adj\n";
			    #$fsadj+=(-1*($f->[4]*(scalar(@{$f->[5]}))));
			    if($f->[4] == 1){
				#substr($encoding,$fmax - $start-1,1,'F');
				substr($gseq,$start-$fmin,1) = '';
				$adj++;
			    }
			    elsif($f->[4] == -1){
				#substr($encoding,$fmax - $start,1,'B');
				#substr($encoding,$fmax- $start,0) = 'B';
				substr($gseq,$start-$fmin+1,0) = 'N';
			    }
			    elsif($f->[4] == 0){
				#substr($encoding,$fmax - $start-1,1,'I');
			    }		
			    die if(($fmax-$start) < 0);
			    #die if(($fmax-$start) >= length($encoding));
			}
		    }
		}
		eval{
		    $newobjs = $newobjs->revcom();
		};
		if($newobjs->seq() =~ /$PMARK_SPACER/){
		    my $sloc = $-[0];
		    print "FOUND PMARK- $sloc ",substr($newobjs->seq(),$sloc,length($PMARK_SPACER)),"\n" if($verbose);
		    #substr($encoding,$sloc,length($PMARK_SPACER),'F'x36);
		    $pmark=1;
		}
		$newobj = new Bio::Seq(-seq=>$gseq);
		$newobj = $newobj->revcom();
		print "GSEQPOST:$gseq\n" if($debug);
		#print "GSEQobj:",$newobj->seq(),"\n";
		
		#print "$seqname $seqlen+$fsadj $fmin,$fmax,$orient PMARK=$pmark ",$newobjs->seq(),"\n",$encoding,"\n";
		print "#Encoding $encoding\n" if($debug);
		#return 0 if(($seqlen+$fsadj)%3!=0);
		#return 0 if(length($encoding) != $newobjs->length());
		#die if(length($encoding) != $newobjs->length());

		# if(0){
# 		    eval{
# 			$newobj = new Bio::Seq::EncodedSeq(-seq=>$newobjs->seq(),
# 							   -encoding=>$encoding,
# 							   -verbose=>0);
# 		    }
# 		    or do{
# 			print "ERROR encoding: ",$@,"\n" if($verbose);
# 			print "$seqname $fmin,$fmax,$orient PMARK=$pmark ",$newobjs->seq(),"\n",$encoding,"\n" if($verbose);
# 			return 0;
# 		    };
# 		}
	    }
	    else{
		$newobj = $seqobj->trunc($fmin+1,$fmax);
		eval{
		    $newobj = $newobj->revcom();
		};
	    }
	    #print "NEW ",$newobj->seq(),"\n";
	    die if($verbose && $newobj->length() > $MAXORFLEN);
	    die if($verbose && $newobj->length() < $MINORFLEN);
	    #Check if valid start codon
	    if(1){#$codon_table->is_start_codon($newobj->subseq(1,3)) && ($codon_table->is_ter_codon($newobj->subseq($newobj->length()-3+1,$newobj->length())))){
		*STDERR=*FOO;
		my $protein_seq_obj;
		eval{
		    if(0 && $fs){
			#print "Using FS\n";
			$newobj->verbose(0);
			$protein_seq_obj = $newobj->cds()->translate(
								     -codontable_id =>11,
								     #-orf=>1,
								     -complete => 1,
								     -throw => 1,
								     -verbose => 0
								     );
		    }
		    else{
			#print "PRETRANS ",$newobj->seq(),"\n";
			$protein_seq_obj = $newobj->translate(-codontable_id =>11,
							      -complete => 1,
							      -throw => 1
							      );
		    }
		}
		or do {
		    print "ERROR translate: ",$@,"\n" if($verbose);
		    print "ERROR translate $seqname $fmin,$fmax,$orient ",$newobj->seq(),"\n" if($verbose);
		    return 0;
		};
		*STDERR=*TEMP;

		return 0 if(!$protein_seq_obj);
		#print $protein_seq_obj->length()," ",$newobj->length()," ",$fsadj,"\n";
		if($protein_seq_obj->length() >0){#== ($newobj->length()+$fsadj)/3){
		    die if($protein_seq_obj->seq() =~ /\*/);
		    return 1+$pmark;
		}
		else{
		    print "#Unexpected sequence length ",$protein_seq_obj->length()," expecting ",($newobj->length()+$fsadj)/3," from ORF $seqname $fmin-$fmax $orient ",$protein_seq_obj->seq(),"\n" if($verbose);
		    for(my $i=0;$i<$protein_seq_obj->length();$i++){
			print $protein_seq_obj->subseq($i+1,$i+1)," ",$newobj->subseq($i*3+1,$i*3+3),"\n";
		    }
		}
	    }
	    else{
		#print "Possible alternative ORF on $seqname $fmin-$fmax,$orient has invalid start:",$newobj->subseq(1,3)," ",$codon_table->is_start_codon($newobj->subseq(1,3))," or stop:",$newobj->subseq($newobj->length()-3+1,$newobj->length())," ",$codon_table->is_ter_codon($newobj->subseq($newobj->length()-3+1,$newobj->length()))," ",$newobj->seq(),"\n" if($verbose);
		return -1;
	    }
	}
    }
    close FOO;
    return 0;
}


#
#Print members and attributes for a cluster
#$query is the longest member of a cluster
#Supported attributes
sub reportCluster{
    my($query,$mappedorgs,$mappedgenes,$unmappedorgs,$unmappedgenes,$feat_attrs,$cluster_attrs,$seq_attrs,$new_orfs) = @_;
    if(scalar(keys %$mappedgenes)>0){
	print $cogfh "COG = $cluster_id, size ",scalar(keys %$mappedgenes), ", connections = 0, perfect = 0;\n";
	print $cogfh "\t$features->{$query}->[5]\n";
	foreach my $organism (sort {$a cmp $b} keys %$mappedorgs){
	    foreach my $gene (sort {$features->{$a}->[1] <=> $features->{$b}->[1]} (keys %{$mappedorgs->{$organism}->{'features'}})){
		if($gene ne $query){
		    print $cogfh "\t$features->{$gene}->[5]\n";
		}
	    }
	}
	my @posscauses = ('CS3','CE3','CS4','CE4','CS0','CE0');
	my $causesstr;
	foreach my $p (@posscauses){
	    if(exists $cluster_attrs->{$p}){
		$causesstr .= "$p;";
		delete $cluster_attrs->{$p};
	    }
	}
	my $classesstr = join(';',sort {$a cmp $b} keys %{$cluster_attrs});
	print ">CLUSTER_$cluster_id num_seqs=",scalar(keys %$mappedorgs)," num_genes=",scalar(keys %$mappedgenes);
	if(exists $mappedgenes->{$query}->{'alignments'}){
	    print " classes=$classesstr query=$query ";
	    if(length($causesstr)>0){
		print " causes=$causesstr ";
	    }
	    print " num_alignments=",scalar(@{$mappedgenes->{$query}->{'alignments'}})," alignments=",join(',',@{$mappedgenes->{$query}->{'alignments'}}) if($debug);
	}
	print "\n";
	print $cfh "CLUSTER_$cluster_id (",scalar(keys %$mappedgenes)," features,",scalar(keys %$mappedorgs)," genomes, classes=$classesstr, query=$query): ";
	print $ctfh "C_$cluster_id\t";
	print $ctfh2 "C_$cluster_id\t";
	$clusters->{$cluster_id}->{'num_feats'} = scalar(keys %$mappedgenes);
	$clusters->{$cluster_id}->{'num_genomes'} = scalar(keys %$mappedorgs);
	$clusters->{$cluster_id}->{'num_alignments'} = scalar(@{$mappedgenes->{$query}->{'alignments'}});
	$clusters->{$cluster_id}->{'classes'} = $classesstr;

	my $qfmin = $features->{$query}->[1];
	my $qfmax = $features->{$query}->[2];
	my $qseqname = $features->{$query}->[0];
	my @mappedfeats;
	my $outtable = [];
	foreach my $organism (sort {$a cmp $b} keys %$mappedorgs){
	    my($start,$end) = &getspan($mappedgenes,keys %{$mappedorgs->{$organism}->{'features'}});
	    my @ogenes = sort {$features->{$a}->[1] <=> $features->{$b}->[1]} (keys %{$mappedorgs->{$organism}->{'features'}});
	    my @ocovs;# = map {sprintf("%.2f",$mappedgenes->{$_}->{'cov'}/$features->{$_}->[3])} (@ogenes); #%coverage over gene length
	    my @oids;#  = map {sprintf("%.2f",$mappedgenes->{$_}->{'pid'}/$mappedgenes->{$_}->{'len'})} (@ogenes); #%id over aligned length
	    my @orients;
	    my @names;
	    my @frames;
	    my @qframes;
	    my @sdist;

	    my $classes;
	    my $longestorf=0;



	    foreach my $gene (@ogenes){
		if(exists $feat_attrs->{$gene}){
		    foreach my $c (sort {$a cmp $b} keys %{$feat_attrs->{$gene}}){
			$classes->{$c}++;
		    }
		}
		$longestorf = ($features->{$gene}->[3] > $longestorf) ? $features->{$gene}->[3] : $longestorf;

		push @ocovs,(sprintf("%.2f",$mappedgenes->{$gene}->{'cov'}/$features->{$gene}->[3])) if($features->{$gene}->[3]>0);
		push @oids,(sprintf("%.2f",$mappedgenes->{$gene}->{'cov'}/$features->{$gene}->[3])) if($features->{$gene}->[3]>0);
		push @orients,"$features->{$gene}->[4]";	
		#Print frame and sdist for inconsistent clusters only
		if(exists $cluster_attrs->{'CS1'} && exists $cluster_attrs->{'CE1'}){
		}
		else{
		    if($mappedgenes->{$gene}->{'relqrysdist'}>0 || $verbose==1){
			push @frames,"altframeqry=$mappedgenes->{$gene}->{'frame'}";
			push @qframes,"frameinqry=$mappedgenes->{$gene}->{'frameinqry'}";
			push @sdist,"sdist=$mappedgenes->{$gene}->{'relqrysdist'}" if(scalar(@ogenes)==1);
		    }
		}
		if(defined $features->{$gene}->[11]){
		    push @names,"product=$features->{$gene}->[11]";
		}
		#push @attrs,"aln_orient=$mappedgenes->{$gene}->{'relorient'}";
		#my $frame;
		#if($features->{$gene}->[4] eq '-'){
		#    $frame=($end%3)*-1;
		#}
		#else{
		#    $frame=$start%3;
		#}
		#push @attrs,"frame=$frame";

	    }

            #Print attributes
	    my @attrs = sort {$a cmp $b} keys %$classes;
	    push @attrs,@frames,@qframes,@sdist;
	    #push @attrs,map {"frame=$_"} @frames;
	    #push @attrs,map {"sdist=$_"} @sdist;


	    #Brief cluster output
	    print $cfh join(',',@ogenes),"($organism) ";
	    my($realorg) = ($organism =~ /([^\.]+)/);
	    $outtable->[$seqindex->{$realorg}]->[0] = join(',',@ogenes);
	    $outtable->[$seqindex->{$realorg}]->[1] = $start;
	    $outtable->[$seqindex->{$realorg}]->[2] = $end;
	    $outtable->[$seqindex->{$realorg}]->[3] = join(',',@orients);

	    #Detailed output
	    print join(',',@ogenes),
	    "\tC$cluster_id",
	    "\t$organism",
	    "\tcov=",join(',',@ocovs),
	    "\tpid=",join(',',@oids),
	    "\tqcov=",sprintf("%.2f",$mappedorgs->{$organism}->{'qcov'}/($qfmax-$qfmin)),
	    "\t$start-$end",
	    "\t",join(',',@orients),
	    "\t",$end-$start,
	    "\t",join(';',@attrs,@names),
	    "\n";

	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'genes'} = \@ogenes;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'cov'} = \@ocovs;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'pid'} = \@oids;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'frame'} = \@frames;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'frameinqry'} = \@qframes;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'sdist'} = \@sdist;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'fmin'} = $start;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'fmax'} = $end;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'len'} = $end-$start;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'orient'} = \@orients;
	    $clusters->{$cluster_id}->{'orgs'}->{$organism}->{'desc'} = join(';',@names);

	}
	if($verbose){
	    foreach my $organism (sort {$a cmp $b} keys %$unmappedorgs){
		if(! exists $unmappedorgs->{$organism}->{'features'}){
		    if(exists $new_orfs->{$organism}){
			print "#ALIGNED_NEWORFS $organism aligned with unannotated matching ORFs\n";
		    }
		    else{
			#No annotated features on this org
			print "#ALIGNED_NOORFS $organism aligned with no matching ORFs, possibly in a gapped region of the alignment\n";
		    }
		}
		else{
		    #print annotated features
		    my($start,$end) = &getspan($unmappedgenes,keys %{$unmappedorgs->{$organism}->{'features'}});
		    my @ogenes = sort {$features->{$a}->[1] <=> $features->{$b}->[1]} (keys %{$unmappedorgs->{$organism}->{'features'}});
		    
		    my @ocovs;
		    my @oids;
		    my @orients;
		    my @names;
		    my $classes;
		    my $longestorf=0;
		    
		    foreach my $gene (@ogenes){
			if(exists $feat_attrs->{$gene}){
			    foreach my $c (sort {$a cmp $b} keys %{$feat_attrs->{$gene}}){
				$classes->{$c}++;
			    }
			}
			$longestorf = ($features->{$gene}->[3] > $longestorf) ? $features->{$gene}->[3] : $longestorf;
			
			push @ocovs,(sprintf("%.2f",$unmappedgenes->{$gene}->{'cov'}/$features->{$gene}->[3])) if($features->{$gene}->[3]>0);
			push @oids,(sprintf("%.2f",$unmappedgenes->{$gene}->{'cov'}/$features->{$gene}->[3])) if($features->{$gene}->[3]>0);
			push @orients,"$features->{$gene}->[4]";		
			
			if(defined $features->{$gene}->[11]){
			    push @names,"product=$features->{$gene}->[11]";
			}
		    }
		    #Print attributes
		    my @attrs = sort {$a cmp $b} keys %$classes;
		    print "#ALIGNED_BELOW_CUTOFFS\t";
		    print join(',',@ogenes),
		    "\tC$cluster_id",
		    "\t$organism",
		    "\tcov=",join(',',@ocovs),
		    "\tpid=",join(',',@oids),
		    "\tqcov=",sprintf("%.2f",$unmappedorgs->{$organism}->{'qcov'}/($qfmax-$qfmin)),
		    "\t$start-$end",
		    "\t",join(',',@orients),
		    "\t",$end-$start,
		    "\t",join(';',@attrs,@names),
		    "\n";
		}
	    }
	}
	    ##
	    #Report ORFs that are conserved and aligned but not annotated
# 	    foreach my $organism (keys %$new_orfs){
# 		my $orfidx=0;
# 		foreach my $alt (@{$new_orfs->{$organism}}){
# 		    die if(exists $mappedorgs->{$organism});
# 		    my($astart,$aend) = ($alt =~ /alt_start=(\d+)-(\d+)/);
# 		    my($len) = ($alt =~ /len:(\d+)/);
# 		    my($orient) = ($alt =~ /orient:([\+\-])/);
# 		    die "Mismatching lengths $len != $aend - $astart" if($len != ($aend-$astart));
# 		    #Check that this ORF is longer than genes that are already annotated on $organism in this region
# 		    my @unmappedlist;
# 		    foreach my $feat_name (keys %$unmappedgenes){
# 			if($feat2organism->{$feat_name} eq $organism){
# 			    push @unmappedlist,[$feat_name,$features->{$feat_name}->[3]];
# 			}
# 		    }
# 		    my @longestunmapped = sort {$b->[1] <=> $a->[1]} @unmappedlist;
# 		    if(scalar (%$unmappedgenes) ==0 || $len > $longestunmapped[0]){
# 			push @mappedfeats,[[[$organism,$astart,$aend,'+',($aend-$astart).'M']],"NEWORF$organism$orfidx",'gene'];
# 			print "NEWORF$organism$orfidx",
# 			"\tWGA$cluster_id",
# 			"\t$organism",
# 			"\tcov=",
# 			"\tpid=",
# 			"\t$astart-$aend",
# 			"\t",$aend-$astart,
# 			"\t",join(';',@{$new_orfs->{$organism}}),
# 			"\n";
# 			$orfidx++;
# 		    }
# 		}
# 	     }

	for(my $i=0;$i<scalar(@sortedallseqs);$i++)
	{
	    print $ctfh $outtable->[$i]->[0];
	    print $ctfh2 "$outtable->[$i]->[0] $outtable->[$i]->[1] $outtable->[$i]->[2] $outtable->[$i]->[3]";
	    if($i!=$#sortedallseqs){
		print $ctfh "\t";
		print $ctfh2 "\t";
	    }
	    else{
		print $ctfh "\n";
		print $ctfh2 "\n";
	    }
	}
	
	if($printalignments){
	    print "#Printing query $query $qseqname,$qfmin,$qfmax\n" if($debug);
	    my $outfh;#=\*STDOUT;
	    open $outfh,"+>$options{'prefix'}cluster_${cluster_id}.aln.out";

	    my @isect = $atree->map($qseqname,$qfmin,$qfmax);
	    #Print all features overlapping the alignment window.
	    #This may include addl features than those in the cluster
	    my $printedfeats = {};
	    foreach my $feat (@isect){
		my $feat_name = $feat->[0];
		$feat_name =~ s/gene\://;
		$printedfeats->{$feat_name}++;
	    }
	    foreach my $feat_name (keys %$printedfeats){
		my $fmin = $features->{$feat_name}->[1];
		my $fmax = $features->{$feat_name}->[2];
		my $seqname = $features->{$feat_name}->[0];
		my $orient = $features->{$feat_name}->[4];
		if(exists $mappedgenes->{$feat_name}){
		    push @mappedfeats,[[[$seqname,$fmin,$fmax,$orient,($fmax-$fmin).'M']],'gene:'.$feat_name,'gene'];
		}
		else{
		    #print "#WARNING Expected gene $feat_name in unmapped list: ".join(',',keys %$unmappedgenes)."\n" if(!exists $unmappedgenes->{$feat_name});
		    if(exists $features->{$feat_name} && $features->{$feat_name}->[3]){
			my $cov = sprintf("c%.1f,i%.1f ",$unmappedgenes->{$feat_name}->{'cov'}/$features->{$feat_name}->[3],
					  $unmappedgenes->{$feat_name}->{'pid'}/$features->{$feat_name}->[3]);
			push @mappedfeats,[[[$seqname,$fmin,$fmax,$orient,($fmax-$fmin).'M']]," $cov *gene:".$feat_name.":$orient",'gene'];
		    }
		}
	    }
	    
	    #Sort all alignments that span query gene
	    my @qryalns;
	    #TODO fix for reverse strand
	    foreach my $align_name (@{$mappedgenes->{$query}->{'alignments'}}){
		die if(!exists $feat2organism->{$query} || length($feat2organism->{$query})==0);
		my $alni = $atree->getAlignedInterval($align_name,$feat2organism->{$query});
		if($alni){
		    print "#QRYALN $align_name $feat2organism->{$query} $alni->[0] $alni->[1] $alni->[2]\n" if($debug);
		    push @qryalns,[$align_name,$alni->[1]];
		}
	    }
	    my $aidx;
	    foreach my $al (sort {$a->[1] <=> $b->[1]} @qryalns){
		my($align_name) = @$al;
		#Check that new range is still within $alignment
		print "#Checking the $qseqname,$qfmin,$qfmax,$align_name is within range\n" if($debug);;
		my @isect = $atree->intersect($qseqname,$qfmin,$qfmax,$align_name);
		my $printfmin;
		my $printfmax;
		foreach my $aln (@isect){
		    if($aln->[1] eq $qseqname && $aln->[0] eq $align_name){
			#print join(',',@$aln),"\n";
			$printfmin = $aln->[2];
			if($aidx==0){
			    if($atree->contains($align_name,$qseqname,$printfmin-20,$printfmin)){
				$printfmin -= 20;
			    }
			}
			$printfmax = $aln->[3];
			if($aidx==(scalar(@qryalns)-1)){
			    if($atree->contains($align_name,$qseqname,$printfmax,$printfmax+20)){
				$printfmax += 20;
			    }
			}
			print "#Resetting print range to $printfmin-$printfmax from $qfmin-$qfmax\n" if($debug);
		    }
		}
		if(defined $printfmin && defined $printfmax){
		    print $outfh "CLUSTER_$cluster_id ALIGNMENT:$align_name\n";
		    my($colstart,$colend) = AlignmentTree::coordstocolumn($atree->{_alignments}->{$align_name}->[0],$qseqname,$printfmin,$printfmax,1);
		    $atree->printAlignment($outfh,$align_name,$colstart,$colend,$db,$qseqname,\@mappedfeats,$htmlout);
		}
		else{
		    die;
		}
		$aidx++;
	    }
	    close $outfh;
	    print `cat $options{'prefix'}cluster_${cluster_id}.aln.out`;

	}
	print "\n";
	print $cfh "\n";
    }
    else{
	#No genes in cluster
	die;
    }
}


#
sub getFeaturesByInterval{
    my($atree,$org,$fmin,$fmax,$orient) = @_;
    my @misects = $atree->intersect($org,$fmin,$fmax,"gene");
    my $feats;
    foreach my $fisectn (@misects){
	#my($fname,$fseq,$fstart,$fend,$fcoverage,$fpid,$forient1,$forient2) = @$fisectn;
	my $feat_name = $fisectn->[0];
	$feat_name =~ s/^gene://;
	$feats->{$feat_name} = [$fisectn->[1],$fisectn->[2],$fisectn->[3],$fisectn->[7],$feat_name,$feat_name,$fisectn->[4],$fisectn->[5]];
    }
    return $feats;
}

#Returns @features,@seqs,@alignments that overlap query
sub getAlignedFeatures{
    my($atree,$seqname,$query,$fmin,$fmax,$type) = @_;
    #Aligned features
    my @nisects;
    #Aligned seqs
    my $seqs = {};
    #Alignments
    my $alignments = {};

    #Parse overlapping genes
    my @isect = $atree->map($seqname,$fmin,$fmax);
    #First screen all overlapping alignments to ensure that they
    #include the query gene
    my $alignments;
    foreach my $r (@isect){
	my $feat_name = $r->[0];
	my $seqname = $r->[1];
	my $align_name = $r->[5];
	#Only consider WGA alignments (alignment name in $align_name) that span query (gene name in $query)
	if($feat_name eq $type.':'.$query){
	    print "#Mapped $query $align_name\n" if($debug);
	    $alignments->{$align_name}++;
	}
    }

    #Capture all seqs in this alignment
    foreach my $align_name (keys %$alignments){
	my $alignedseqs  = $atree->{_alignments}->{$align_name}->[0];
	foreach my $seq (@$alignedseqs){
	    die if(ref $seq->[0]);
	    $seqs->{$seq->[0]}++;
	    #TODO capture stats on alignment
	    #$seqs->{$seq->[0]}->{'len'}+=;
	    #$seqs->{$seq->[0]}->{'pid'}+=;
	    #$seqs->{$seq->[0]}->{'cov'}+=;
	}
    }

    #Transform feat_name, stripping leading type:
    my @nisect;
    foreach my $r (@isect) {
	$r->[0] =~ s/$type\://;
	if(! exists $features->{$r->[0]} && !defined $featlist){
	    print STDERR "Unknown feature $r->[0]\n";
	}
	print "#SAM$r->[0] ",(exists $features->{$r->[0]}),"\n" if($debug);
	push @nisect,$r if(exists $features->{$r->[0]});
    }
    
    return (\@nisect,$seqs,$alignments);    
}

sub findSingletons{
    my($atree,$mapped,$unmapped,$subsumed,$datree) = @_;
    my $singletons = {};
    my $dups = {};
    foreach my $feat_name (keys %$features){
	my $fmin = $features->{$feat_name}->[1];
	my $fmax = $features->{$feat_name}->[2];
	my $seqname = $features->{$feat_name}->[0];
	if(! exists $mapped->{$feat_name}){
	    die if(exists $mapped->{$feat_name});
	    my $classes = &annotateSingletons($atree,$features->{$feat_name}->[0],$feat_name,$fmin,$fmax);

	    my $nisect;
	    my $allseqs;
	    my $goodalignments;
	    my $dupfeats;
	    if(defined $datree){
		#print "Querying $seqname,$feat_name,$fmin,$fmax,'gene'\n";
		($nisect,$allseqs,$goodalignments) = &getAlignedFeatures($datree,$seqname,$feat_name,$fmin,$fmax,'gene');
		#print scalar(@$nisect),"\n";
		foreach my $r (
			       sort { $features->{$b->[0]}->[3] <=> $features->{$a->[0]}->[3] } #sort on feature length
			       @$nisect){
		    my $dfeat_name = $r->[0];
		    my $seqname = $r->[1];
		    my $align_name = $r->[5];
		    #Check if we want to consider this alignment
		    if(exists $goodalignments->{$align_name}){
			$dfeat_name =~ s/gene\://;
			if(!exists $features->{$dfeat_name}){
			    print STDERR "#Bad feature found $dfeat_name. Not in input file. Skipping\n";
			    next;
			}
			if($dfeat_name ne $feat_name){
			    #print "Saving dup $dfeat_name , $feat_name\n";
			    $dupfeats->{$feat_name}->{$dfeat_name}->{'cov'} += $r->[7];
			    $dupfeats->{$feat_name}->{$dfeat_name}->{'pid'} += $r->[8];
			    $dupfeats->{$feat_name}->{$dfeat_name}->{'len'} += ($r->[3]-$r->[2]);
			}
		    }
		}
	    }

	    if(exists $unmapped->{$feat_name}){
		my $query=$feat_name;
		my($mappedorgs,$mappedgenes,$unmappedorgs,$unmappedgenes) = &buildCluster($atree,$query);
		my($feat_attrs,$cluster_attrs,$codons) = &annotateCluster($atree,$mappedgenes,$mappedorgs);
		my $new_orfs = &findnewORFs($db,$atree,$mappedorgs,$mappedgenes,$codons);
		if(scalar(keys %$new_orfs)){
		    my $seq_attrs = {};	 
		    &reportCluster($query,$mappedorgs,$mappedgenes,$unmappedorgs,$unmappedgenes,$feat_attrs,$cluster_attrs,$seq_attrs,$new_orfs);
		    $cluster_id++;
		}
		my $featlen = $fmax-$fmin;
		my $mappedlen = $unmapped->{$feat_name}->{'len'};
		if($featlen <= 0){
		    print STDERR "#Bad featlen for feature $feat_name $fmax-$fmin\n";
		    $featlen=1;
		}
		if($mappedlen <= 0){
		    print STDERR "#Bad coverage for feature $feat_name Coverage:$unmapped->{$feat_name}->{'len'}\n";
		    $mappedlen=1;
		}
		my ($seqname,$fmin,$fmax,$len,$orient) = @{$features->{$feat_name}};
		print "#SINGLETON $feat_name len:$features->{$feat_name}->[3]\tbest_cluster:C$unmapped->{$feat_name}->{'WGA_cluster'}\tcov:";
		#printf("%.2f",$unmapped->{$feat_name}->{'cov'}/$featlen);
		printf("%.2f",$unmapped->{$feat_name}->{'cov'});
		print " pid:";
		#printf("%.2f",$unmapped->{$_}->{'pid'}/$mappedlen);
		printf("%.2f",$unmapped->{$feat_name}->{'pid'});
		printf(" lenbp:%f ",$mappedlen);
		join(' ',@$classes);
		if(defined $features->{$feat_name}->[11]){
		    print " product=$features->{$feat_name}->[11]";
		}
		if(scalar(keys %{$dupfeats->{$feat_name}}) > 0){
		    foreach my $dfeat_name (sort {$dupfeats->{$feat_name}->{$b}->{'pid'} <=> $dupfeats->{$feat_name}->{$a}->{'pid'}} keys %{$dupfeats->{$feat_name}}){
			print " #DUP matches:$dfeat_name(pid:";
			    printf("%.2f",($dupfeats->{$feat_name}->{$dfeat_name}->{'pid'}/$features->{$feat_name}->[3]));
			print ",cov:";
			printf("%.2f",($dupfeats->{$feat_name}->{$dfeat_name}->{'cov'}/$features->{$feat_name}->[3]));
			print ") ";
		    }
		}
		print "\n";
	    }
	    else{
		if(exists $subsumed->{$feat_name}){
		    print "#DELETED $feat_name\n";
		}
		else{
		    print "#SINGLETON $feat_name len:$features->{$feat_name}->[3] ",join(' ',@$classes);
		    if(defined $features->{$feat_name}->[11]){
			print " product=$features->{$feat_name}->[11]";
		    }
		    if(scalar(keys %{$dupfeats->{$feat_name}}) > 0){
			foreach my $dfeat_name (keys %{$dupfeats->{$feat_name}}){
			    print " #DUP matches:$dfeat_name(pid:";
				printf("%.2f",($dupfeats->{$feat_name}->{$dfeat_name}->{'pid'}/$features->{$feat_name}->[3]));
			    print ",cov:";
			    printf("%.2f",($dupfeats->{$feat_name}->{$dfeat_name}->{'cov'}/$features->{$feat_name}->[3]));
			    print ") ";
			}
			$dups->{$feat_name}++;
		    }
		    else{
			$singletons->{$feat_name}++;
		    }
		    print "\n";
		}
	    }
	}
	else{
	    #Mapped ORF, not a singleton
	}
    }

    return ($singletons,$dups);
}


###############################
#General utility funcs
sub getspan{
    my($features) = shift;
    my @coords;
    foreach my $gene (@_){
	die if(! exists $features->{$gene});
	die if(! exists $features->{$gene});
	push @coords,$features->{$gene}->{'fmin'},$features->{$gene}->{'fmax'};
    }
    my @sortedcoords = sort {$a <=> $b} @coords;
    return ($sortedcoords[0],$sortedcoords[$#coords]);
}


sub findCoords{
    my($atree,$seqname,$startcodon,$stopcodon) = @_;
    
    #$codon is a tuple of alignment,aligned_column
    my($startcol,$aln_s) = split(/$CODON_DELIM_REGEX/,$startcodon);
    #find corresponding stop
    my($stopcol,$aln_e) = split(/$CODON_DELIM_REGEX/,$stopcodon);
    
    my $si = &getAlignment($atree,$aln_s,$seqname);
    my $ei = &getAlignment($atree,$aln_e,$seqname);
    my $start_s;
    my $start_e;
    my $stop_s;
    my $stop_e;
    if($si){
	($start_s,$start_e) = AlignmentTree::columntocoords($si,$startcol,$startcol+2);
	if($start_s == $start_e){
	    #aligned to a gap
	    $start_s = undef;
	    $start_e = undef;
	}
	if($ei){
	    ($stop_s,$stop_e) = AlignmentTree::columntocoords($ei,$stopcol,$stopcol+2);
	    if($stop_s == $stop_e){
		#aligned to a gap
		$stop_s = undef;
		$stop_e = undef;
	    }
	}
	else{
	    print "Can't find alignment $aln_s on $seqname from $startcodon\n" if($debug);
	    return undef;
	}
    }
    else{
	print "Can't find alignment $aln_s on $seqname from $startcodon\n" if($debug);
	return undef;
    }
    if($start_s<$stop_s){
	#forward strand 5'start -----> 3'stop
	return ($start_s,$stop_e,'+');
    }
    else{
	#reverse strand
	#3'stop <-- 5'start
	return ($stop_s,$start_e,'-');
    }
    return undef;
}

#Returns aligned location of start and stop codons
#If annotation is not a valid start or stop codons returns -1
#If codon is not aligned returns undef
sub findCodons{
    my($atree,$seqname,$fmin,$fmax,$orient,$fname) = @_;
    #my($name,$seq,$start,$end,$coverage,$qpid) = @$aln;
    my $codon_table = Bio::Tools::CodonTable->new(-id=>11);
    #my $seqobj = $db->get_Seq_by_id($seqname);
    my $seqobj = $db->{$seqname};
    if(!$seqobj){
	print "Can't find seqname: $seqname\n";
	return;
    }
    my $startcodon=undef;
    my $stopcodon=undef;
    my $is_partial_start=0;
    my $is_partial_stop=0;
    my $is_bad_start=0;
    my $is_bad_stop=0;
    my $aln_orient=undef;
    if($orient eq '+'){
	if($fmin+1<=0){print STDERR "Bad start parameter $fmin+1<=0 $seqname,$fmin,$fmax,$orient,$fname\n";return} 
	if($fmax-3+1<=0){print STDERR "Bad end parameter $fmax-3+1<=0 $seqname,$fmin,$fmax,$orient,$fname\n";return}
	if($fmin+3>$seqobj->length()){print STDERR "Bad start parameter $fmin+1<=0 $seqname,$fmin,$fmax,$orient,$fname\n";return};
	if($fmax>$seqobj->length()){print STDERR "Bad end parameter $fmax-3+1<=0 $seqname,$fmin,$fmax,$orient,$fname\n";return};
	
	if(!$codon_table->is_start_codon($seqobj->subseq($fmin+1,$fmin+3))){ #bioperl is 1-base coordinates
	    print "#Bad start codon $fname,$seqname,$fmin,$fmax,$orient codon $fmin+1,$fmin+2+1 ",$seqobj->subseq($fmin+1,$fmin+3)," aln_orient:$aln_orient\n" if($verbose || $debug);
	    $startcodon = &getAlignedCols($atree,$seqname,$fmin,$fmin+3);
	    $is_bad_start=1;
	} 
	else{
	    #Find start codon + strand
	    $startcodon = &getAlignedCols($atree,$seqname,$fmin,$fmin+3);
	}
	
	if(!$codon_table->is_ter_codon($seqobj->subseq($fmax-3+1,$fmax))){
	    print "#Bad stop $fname,$seqname,$fmin,$fmax,$orient codon $fmax-3+1,$fmax ",$seqobj->subseq($fmax-3+1,$fmax)," aln_orient:$aln_orient\n" if($verbose || $debug);
	    $stopcodon = &getAlignedCols($atree,$seqname,$fmax-3,$fmax);
	    $is_bad_stop=1;
	}
	else{
	    #Find stop codon - strand
	    $stopcodon = &getAlignedCols($atree,$seqname,$fmax-3,$fmax);
	}
	#Check if in pmark spacer adjacent to contig boundary
	if($fmin-length($PMARK_SPACER)>0 && $fmin+length($PMARK_SPACER) <= $seqobj->length()){
	    my $startregion = $seqobj->subseq($fmin-length($PMARK_SPACER),$fmin+length($PMARK_SPACER));
	    if($startregion =~ /$PMARK_SPACER/){
		$is_partial_start=1;
	    }
	}
	if($fmax-length($PMARK_SPACER)>0 && $fmax+length($PMARK_SPACER) <= $seqobj->length()){
	    my $stopregion = $seqobj->subseq($fmax-length($PMARK_SPACER),$fmax+length($PMARK_SPACER));
	    if($stopregion =~ /$PMARK_SPACER/){
		$is_partial_stop=1;
	    }
	}
	
    }
    else{
	die "Bad orient $orient" if($orient ne '-');
	print STDERR "Bad start parameter $fmin+1<=0 $seqname,$fmin,$fmax,$orient,$fname" if($fmax-3+1<=0);
	print STDERR "Bad end parameter $fmax-3+1<=0 $seqname,$fmin,$fmax,$orient,$fname" if($fmin+1<=0);
	print STDERR "Bad start parameter $fmax>$seqobj->length() $seqname,$fmin,$fmax,$orient,$fname" if($fmax>$seqobj->length());
	print STDERR "Bad end parameter $fmin+3>$seqobj->length() $seqname,$fmin,$fmax,$orient,$fname" if($fmin+3>$seqobj->length());
	eval{
	    if(!$codon_table->is_start_codon(revcom($seqobj->subseq($fmax-3+1,$fmax))->seq())){
		print "#Bad start codon $fname,$seqname,$fmin,$fmax,$orient codon $fmax-3+1,$fmax ",revcom($seqobj->subseq($fmax-3+1,$fmax))->seq()," aln_orient:$aln_orient\n" if($verbose || $debug);
		$startcodon = &getAlignedCols($atree,$seqname,$fmax-3,$fmax);
		$is_bad_start=0;
	    } 
	    else{
		#Find start codon on - strand
		$startcodon = &getAlignedCols($atree,$seqname,$fmax-3,$fmax);
	    }
	    if(!$codon_table->is_ter_codon(revcom($seqobj->subseq($fmin+1,$fmin+3))->seq())){
		print "#Bad stop codon $fname,$seqname,$fmin,$fmax,$orient codon $fmin+1,$fmin+3 ",revcom($seqobj->subseq($fmin+1,$fmin+3))->seq()," aln_orient:$aln_orient\n" if($verbose || $debug);
		$stopcodon = &getAlignedCols($atree,$seqname,$fmin,$fmin+3);
		$is_bad_stop=1;
	    } 
	    else{
		#Find stop codon on - strand
		$stopcodon = &getAlignedCols($atree,$seqname,$fmin,$fmin+3);
	    }
	
	    #Check if in pmark spacer adjacent to contig boundary
	} or do{
	    warn $@ if($verbose);
	    print STDERR "ERROR invalid start,stop codons or invalid translation. $seqname,$fmin,$fmax,$orient,$fname\n";
	    return 0;
	};
	if($fmax-length($PMARK_SPACER)>0 && $fmax+length($PMARK_SPACER) <= $seqobj->length()){
	    my $startregion = $seqobj->subseq($fmax-length($PMARK_SPACER),$fmax+length($PMARK_SPACER));
	    if($startregion =~ /$PMARK_SPACER/){
		$is_partial_start=1;
	    }
	}
	if($fmin-length($PMARK_SPACER)>0 && $fmin+length($PMARK_SPACER) <= $seqobj->length()){
	    my $stopregion = $seqobj->subseq($fmin-length($PMARK_SPACER),$fmin+length($PMARK_SPACER));
	    if($stopregion =~ /$PMARK_SPACER/){
		$is_partial_stop=1;
	    }
	}
    }
    return ($startcodon,$stopcodon,$is_partial_start,$is_partial_stop,$is_bad_start,$is_bad_stop);
}


sub getAlignment{
    my($atree,$align_name,$seqname) = @_;
    my $alignment = $atree->{_alignments}->{$align_name}->[0];
    foreach my $i (@$alignment){
	if($i->[0] eq $seqname){
	    return $i;
	}
    }
    print "#Can't find $seqname on alignment $align_name\n" if($debug);
    return undef;
}

#Look for indels in alignment columns [$codon-$offset,$codon+2]
#Refseq is optional, otherwise uses most frequently occuring allele as reference
#Returns
#[coord,refchar,qrychar,column,frame]
sub reportVariants{
    my($atree,$db,$aln,$seq,$startcol,$endcol,$refseq) = @_;
    my $skipgapcheck=0;
    my $GAPWINDOW=10;
    die if($endcol<$startcol);
    print "#Analyzing codon position $startcol in alignment $aln seq $seq \n" if($debug);

    print "#Retrieving alignment matrix for $startcol-$endcol for alignment $aln \n" if($debug);

    my ($mmatrix,$seqmatrix,$names) = $atree->getAlignmentMatrix($aln,$startcol,$endcol,$db,$refseq,$seq);
    print "#Expecting width ",($endcol-$startcol+1)," row count ",scalar(@$mmatrix)," ",scalar(@$names),"\n" if($debug);
    #List of columns with variants
    my $results = {};
    my @edits;

    my $qryidx;
    #For optional reference seq
    my $refidx=-1;

    my $width;
    for(my $i=0;$i<@$mmatrix;$i++){
	if($names->[$i] eq $seq){
	    $qryidx = $i;
	}
	if(defined $refseq && $names->[$i] eq $refseq){
	    $refidx = $i;
	}
    }
    #Matrix cols start at 0
    for(my $j=0;$j<($endcol-$startcol+1);$j++){
	if(defined $refseq){
	    if(uc(substr($mmatrix->[$refidx],$j,1)) ne uc(substr($mmatrix->[$qryidx],$j,1))){
		#print uc(substr($mmatrix->[$refidx],$j,1))," ", uc(substr($mmatrix->[$qryidx],$j,1)),"\n";
		$results->{$j}++;
	    }
	}
	else{
	    for(my $i=0;$i<@$mmatrix;$i++){
		if(substr($mmatrix->[$i],$j,1) ne '.'){
		    if($skipgapcheck || substr($mmatrix->[$i],$j,$GAPWINDOW) =~ /\./ ){ #gap < GAPWINDOW
			#column $i has multiple characters, gaps or mutations
			print "#MUT $i $j ",substr($mmatrix->[$i],$j,1)," $names->[$i] $seq\n" if($debug);
			$results->{$j}++;
		    }
		}
	    }
	}
    }
    my $alni = &getAlignment($atree,$aln,$seq);
    foreach my $r (sort {$a cmp $b} keys %$results){
	my $reloffset = $startcol+$r;
	my $freqchar = {};
	my $refchar;
	my $qrychar;
	if(defined $refseq){
	    $qrychar = substr($mmatrix->[$qryidx],$r,1);
	    $refchar = substr($mmatrix->[$refidx],$r,1);
	}
	else{
	    for(my $i=0;$i<@$mmatrix;$i++){
		
		my $char;
		#TODO this is slow, improve perf
		if(substr($mmatrix->[$i],$r,1) eq '-'){
		    #gap
		    $char = substr($mmatrix->[$i],$r,1);
		    #die "Unexpected char $i $r $seqmatrix->[$i]->[$r] $mmatrix->[$i]->[$r]" if(defined $seqmatrix->[$i]->[$r]);
		}
		else{
		    #retrieve base
		    $char = substr($seqmatrix->[$i],$r,1);
		}
		die "Bad char '$char'" if(length($char)!=1);
		$freqchar->{$char}++;
	    }
	}



	my $fstype = 0;

	if(defined $refseq){
	    if(uc($refchar) ne uc($qrychar)){
		if($refchar eq '-'){
		    $fstype=1;
		}
		elsif($qrychar eq '-'){
		    $fstype=-1;
		}
		else{
		    $fstype=0;
		}
		#Ignore point mutations for now
		if($fstype!=0){
		    my($fsstart,$fsend) = AlignmentTree::columntocoords($alni,$reloffset,$reloffset);		
		    if(1){#$fsstart != $fsend){
			print "#ALT col:$reloffset coord:$fsstart-$fsend base:$refchar freq:$freqchar->{$refchar} $seq:$qrychar $freqchar->{$qrychar} fstype:$fstype\n" if($debug);
			push @edits,[$fsstart,$refchar,$qrychar,$reloffset,$fstype];
			#if(scalar(@edits)>$FS_THRESHOLD){
			#return \@edits;
			#}
		    }
		}
	    }
	}
	else{
	    die;
	    #report most frequent character
	    my @sortedchars = sort {$b <=> $a} (keys %$freqchar);
	    #retrieve coordinate on $seq for reloffset
	    foreach my $base (@sortedchars){
		if(uc($base) ne uc($qrychar)
		   #&& $freqchar->{$base}>=$freqchar->{$qrychar}	    #only consider bases that occur more frequently than 
		   #&& $freqchar->{$base}>=scalar(@$mmatrix)/2){  	    #optionally also in majority of sequences
		   ){
		    if($base eq '-'){
			$fstype=1;
		    }
		    elsif($qrychar eq '-'){
			$fstype=-1;
		    }
		    else{
			$fstype=0;
		    }
		    my($fsstart,$fsend) = AlignmentTree::columntocoords($alni,$reloffset,$reloffset);		    
		    print "#ALT col:$reloffset coord:$fsstart-$fsend base:$base freq:$freqchar->{$base} $seq:$qrychar $freqchar->{$qrychar} fstype:$fstype\n" if($debug);
		    push @edits,[$fsstart,$base,$qrychar,$reloffset,$fstype];
		}
		else{
		    #last;#can shortcircuit, only consider more frequent bases
		}
	    }
	}
    }
    return \@edits;
}

#Returns the overlapping alignment and start-end column for a sequence range
#Inputs
#getAlignedCols(seq,fmin,fmax)
#Returns [start_colnum,alignment_obj,end_colnum,matching_bits]
sub getAlignedCols{
    my($atree,$seqname,$fmin,$fmax) = @_;
    my $ret;
    my @alignments = $atree->intersect($seqname,$fmin,$fmax,$aligntoken);
    my $found=0;
    foreach my $aln (@alignments){
	if($seqname eq $aln->[1]){
	    my $align_name = $aln->[0];
	    my $align_start = $aln->[2];
	    my $align_end = $aln->[3];
	    die "Bad alignment name $align_name" if(!exists $atree->{_alignments}->{$align_name});
	    die "Mis-mathed orient $aln->[6] ne $aln->[7]" if($aln->[7] ne $aln->[6]);
	    my $alni = $atree->{_alignments}->{$align_name}->[0];
	    if($align_start == $fmin && $fmax == $align_end){
		if($found){
		    print "#WARNING Overlapping aligned region found for $seqname,$fmin,$fmax. $align_name and $ret->[1]\n" if($debug);
		}
		my @res= AlignmentTree::coordstocolumn($alni,$seqname,$fmin,$fmax,1);
		$ret = [$res[0],$align_name,$res[1],$res[2]];
		$found=1;
	    }
	}
    }
    return $ret;
}

################
#DEPRECATED CODE
##############################
#
#callORF()
#Attempts to call an ORF using start codon specified by [start-end]
#Start,end should be codon coordinates relative to the + strand. start<end

#Will attempt to call an ORF on one strand.
#Leading strand 5'->3' increasing coordinates [start-firstStop] 
#Lagging strand 5'->3' decreasing coordinates [end-firstStop]

#Will only call ORF if start,end,orient corresponds to an acutal start
#codon, specified by the configurable codon table

#fsedits is a array reference of signed locations of the frameshift relative to the sequence start
#eg. +10 is a forward frameshift 10 bp downstream from translation start
#    -9 is a backward frameshift 9 bp downstream from translation start
sub callORF{
    my($seqobj,$codon_start,$codon_end,$orient,$fs) = @_;
    die "Bad start codon $seqobj:$codon_start-$codon_end $orient" if($codon_end < $codon_start || $codon_end - $codon_start != 3);
    my $codon_table = Bio::Tools::CodonTable->new(-id=>11);
    if($seqobj){
	if($orient eq '+'){
            my $seqlen = ($seqobj->length()>$MAXORFLEN) ? $codon_start+$MAXORFLEN : $seqobj->length(); 
	    my $newobjs = $seqobj->trunc($codon_start+1,$seqlen);
	    my $encoding = 'C'x$newobjs->length();

	    foreach my $fs_loc (@$fs){
		if(defined $fs_loc){
		    if($fs_loc>0){
			print "Encoding a forward frameshift at $fs_loc in ORF of length ",$newobjs->length(),"\n";
			#a forward frameshift
			#substr($encoding,$fs_loc,1) = 'F';
			substr($encoding,$fs_loc,1,'F');
		    }
		    else{
			#a backward frameshift
			print "Encoding a reverse frameshift at $fs_loc in ORF of length ",$newobjs->length(),"\n";
			#substr($encoding,($fs_loc*-1),1) = 'B';
			substr($encoding,($fs_loc*-1),1,'B');
		    }
		}
	    }
	    die if(length($encoding)!=$newobjs->length());
	    my $newobj = new Bio::Seq::EncodedSeq(-seq=>$newobjs->seq(),
						  -encoding=>$encoding);


	    #Check if valid start codon
	    if($codon_table->is_start_codon($newobj->subseq(1,3))){
		my $protein_seq_obj = $newobj->translate(-orf => 1,
							 -codontable_id =>11);
		return ($protein_seq_obj->seq(),$orient);
	    }
	    else{
		print "#callORF trying '-' $seqobj,$codon_start,$codon_end,$orient Bad start codon ",$newobj->subseq(1,3) if($debug);;
		my $seqlen = ($codon_end>$MAXORFLEN) ? $codon_end-$MAXORFLEN : 1;
		my $newobj = $seqobj->trunc($seqlen,$codon_end);
		eval{
		    $newobj = $newobj->revcom();
		};
		#print " REV:",$codon_table->is_start_codon($newobj->subseq(1,3))," ",$newobj->subseq(1,3),"\n";
		if($codon_table->is_start_codon($newobj->subseq(1,3))){
		    my $protein_seq_obj = $newobj->translate(-orf => 1,
							     -codontable_id =>11);
		    
		    return ($protein_seq_obj->seq(),'-');
		}
		else{
		    print "#WARNING: Skipping callORF $seqobj,$codon_start,$codon_end,$orient. '",$newobj->subseq(1,3),"' is not a valid start codon\n" if($debug);
		}
	    }		
	}
	else{
	    die if($orient ne '-');
            my $seqlen = ($codon_end>$MAXORFLEN) ? $codon_end-$MAXORFLEN : 1;
	    my $newobj = $seqobj->trunc($seqlen,$codon_end);
	    eval{
		$newobj = $newobj->revcom();
	    };
	    #Check if valid start codon
	    if($codon_table->is_start_codon($newobj->subseq(1,3))){
		my $protein_seq_obj = $newobj->translate(-orf => 1,
							 -codontable_id =>11);
		
		return ($protein_seq_obj->seq(),$orient);
	    }
	    else{
		print "#callORF trying '+' $seqobj,$codon_start,$codon_end,$orient Bad start codon ",$newobj->subseq(1,3) if($debug);;
		my $seqlen = ($seqobj->length()>$MAXORFLEN) ? $codon_start+$MAXORFLEN : $seqobj->length(); 
		my $newobj = $seqobj->trunc($codon_start+1,$seqlen);
		if($codon_table->is_start_codon($newobj->subseq(1,3))){
		    my $protein_seq_obj = $newobj->translate(-orf => 1,
							     -codontable_id =>11);
		    
		    return ($protein_seq_obj->seq(),'+');
		}
		else{
		    print "WARNING: Skipping callORF $seqobj,$codon_start,$codon_end,$orient. '",$newobj->subseq(1,3),"' is not a valid start codon\n";
		}
	    }
	}
	    
    }
    else{
	print "#ERROR invalid seq obj $seqobj\n" if($debug);;
    }
    return undef;
}

#Print alternative start sites
#
#Method:
#Report aligned but un-annotated start codons
#Reports
#(1) alternative start location, frequency annotated in the alignment
#(2) resulting ORF, len
sub checkStarts{
    my ($db,$codons,$seqs,$seq_attrs) = @_;

    my $altorfs;
    #Save list of all start codons $codon->$freq
    my $starts = {};
    my $stops = {};
    die if(!exists $codons->{'starts'});
    foreach my $seqname (keys %{$codons->{'starts'}}){
	foreach my $codon (keys %{$codons->{'starts'}->{$seqname}}){
	    $starts->{$codon} += $codons->{'starts'}->{$seqname}->{$codon};
	}
    }
    die if(!exists $codons->{'stops'});
    foreach my $seqname (keys %{$codons->{'stops'}}){
	foreach my $codon (keys %{$codons->{'stops'}->{$seqname}}){
	    $stops->{$codon} += $codons->{'stops'}->{$seqname}->{$codon};
	}
    }

    foreach my $seqname (@{$seqs}){
	#Consider all codons that are not currently annotated on this sequence
	foreach my $codon (keys %$starts){
	    if(! exists $codons->{'starts'}->{$seqname}->{$codon}){ #start codon is not annotated on $seqname
		print "#CODON $codon not annotated on $seqname\n" if($debug);;
		#check if $codon is aligned
		#$codon is a tuple of alignment,aligned_column
		my($col,$aln) = split(/$CODON_DELIM_REGEX/,$codon);
		my $gapped=1; #isgapped
		#check is $col,$col+3 is gapped, return start coordinate on the genome
		my $i = &getAlignment($atree,$aln,$seqname);
		if($i){
		    die "Cannot find alignment $aln that contains $seqname" if(!$i);
		    #Obtain coordinates of the putative start codon
		    my($start,$end) = AlignmentTree::columntocoords($i,$col,$col+2);
		    $gapped = (abs($end-$start) == 3) ? 0 : 1;
		    if(!$gapped){
			#codon is aligned, attempt to call ORF
			#save and report it. save frequency
			if($db){
			    my $orient = $i->[3];
			    print "#Looking for ORF $start,$end,$orient on $seqname\n" if($debug);
			    #my $seqobj = $db->get_Seq_by_id($seqname);
			    my $seqobj = $db->{$seqname};
			    if($seqobj){
				die "Can't find sequence $seqname obj:$seqobj" if(!defined $seqobj);
				my ($neworf,$callorient) = &callORF($seqobj,$start,$end,$orient);
				if(length($neworf)>$MINORFLEN){
				    print "#Calling ORF on strand $callorient start coord = $start\n" if($debug);;
				    #$codons->{'alt_starts'}->{$seqname}->{$codon}->{'freq'} = $starts->{$codon};
				    #$codons->{'alt_starts'}->{$seqname}->{$codon}->{'neworf'} = $neworf;
				    #$codons->{'alt_starts'}->{$seqname}->{$codon}->{'orient'} = $callorient;
				    my $fmin;
				    my $fmax;
				    if($callorient eq '+'){
					$fmin=$start;
					$fmax=$start+(length($neworf)*3);
				    }
				    else{
					$fmin=$end-(length($neworf)*3);
					$fmax=$end;
				    }
				    if(!$fmin || $fmin<0){
					print STDERR "Bad ORF call on $seqname $start,$end converted to $fmin,$fmax\n";
					next;
				    }
				    #$codons->{'alt_starts'}->{$seqname}->{$codon}->{'start'} = $fmin;
				    #$codons->{'alt_starts'}->{$seqname}->{$codon}->{'end'} = $fmax;
				    my($strc,$stpc) = &findCodons($atree,
								  $seqname,
								  $fmin,
								  $fmax,
								  $callorient);
				    #if($callorient eq '-'){
				#	($strc,$stpc) = ($stpc,$strc);
					
				#    }
				    my $startcodon;
				    my $stopcodon;
				    if(ref $strc){
					my($mcol,$align_name) = (@$strc);
					$startcodon = $mcol.$CODON_DELIM.$align_name;
					#die "Can't find start $mcol,$align_name $callorient,$orient from $seqname $codon" if(!exists $starts->{$startcodon});
					#if(!exists $starts->{$startcodon}){
					#    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'startfreq'} = 0;
					#}
					#else{
					#    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'startfreq'} = $starts->{$startcodon};
					#}
					#$codons->{'alt_starts'}->{$seqname}->{$codon}->{'startcol'} = $mcol;
					#$codons->{'alt_starts'}->{$seqname}->{$codon}->{'startcodon'} = $startcodon;
				    }
				    if(ref $stpc){
					my($mcol,$align_name) = (@$stpc);
					$stopcodon = $mcol.$CODON_DELIM.$align_name;
					#die "Can't find stop $mcol,$align_name $callorient,$orient from $seqname $codon" if(!exists $stops->{$stopcodon});
					#if(!exists $stops->{$stopcodon}){
					#    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'stopfreq'} = 0;
					#}
					#else{
					#    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'stopfreq'} = $stops->{$stopcodon};					
					
				        #}
					#$codons->{'alt_starts'}->{$seqname}->{$codon}->{'stopcol'} = $mcol;
					#$codons->{'alt_starts'}->{$seqname}->{$codon}->{'stopcodon'} = $stopcodon;
				    }
				    #Save start,stop pair
				    if($startcodon && $stopcodon){
					#$codons->{'pairs'}->{$startcodon.':'.$stopcodon}->{'gfreq'}++;
					#$codons->{'pairs'}->{$startcodon.':'.$stopcodon}->{'length'} += ($fmax-$fmin);
					#$codons->{'pairs'}->{$startcodon.':'.$stopcodon}->{'orgs'}->{$seqname} = [$fmin,$fmax,0];
					push @$altorfs,[$seqname,$fmin,$fmax,$callorient,$startcodon,$stopcodon];
				    }
				}
				else{
				    print "Skipping short ORF ",length($neworf)," <$MINORFLEN $start,$end,$orient\n" if($debug);
				}
			    }
			    else{
				print "#WARNING. Sequence $seqname not found in FASTA file. Skipping calling new ORFs.\n";
			    }
			}
			else{
			    print "#WARNING. No FASTA file, cannot call new ORFs\n";
			}
		    }
		    else{
			print "#alignment to codon contains gaps $start,$end\n" if($debug);;
		    }
		}
		else{
		    print "#can't find alignment $aln $seqname\n" if($debug);;

		}
	    }
	}
    }
    return $altorfs;
}


sub findNearestNeighbor{
    my($atree,$seqname,$mappedseqs,$start,$end) = @_;
    my @res = $atree->map($seqname,$start,$end);
    my @sres = sort {$b->[8] <=> $a->[8]} @res;
    foreach my $s (@sres){
	if($s->[1] ne $seqname && exists $mappedseqs->{$s->[1]}){
	    return $s->[1];
	}
    }
}

sub reportFrameShifts{
    my($atree,$db,$seqname,$nearestseq,$startcodon,$stopcodon) = @_;
    #$codon is a tuple of alignment,aligned_column
    my($startcol,$aln_s) = split(/$CODON_DELIM_REGEX/,$startcodon);
    #find corresponding stop
    my($stopcol,$aln_e) = split(/$CODON_DELIM_REGEX/,$stopcodon);
    my $si = &getAlignment($atree,$aln_s,$seqname);
    my $ei = &getAlignment($atree,$aln_e,$seqname);
    my($startcoord) = AlignmentTree::columntocoords($si,$startcol,$startcol);
    my($stopcoord) = AlignmentTree::columntocoords($ei,$stopcol,$stopcol);
    #Make sure we have not traversed a rearrangement
    if(abs($stopcoord-$startcoord)<$MAXORFLEN){
	
	my $fsvars = [];
	my $netfs = 0;
		
	#TODO, relax to allow multiple spanning alignments
	print "#Looking for frameshifts in $startcodon,$stopcodon $aln_s $aln_e $startcoord $stopcoord\n" if($debug);
	my @sortedproj;
	if($startcoord < $stopcoord){
	    my @proj = $atree->intersect($seqname,$startcoord,$stopcoord,'WGA');
	    @sortedproj = sort {$a->[2] <=> $b->[2]} @proj;
	}
	else{
	    my @proj = $atree->intersect($seqname,$stopcoord,$startcoord,'WGA');
	    @sortedproj = sort {$b->[3] <=> $a->[3]} @proj;
	}
	print "#Found ",scalar(@sortedproj)," alignments\n" if($debug);
	foreach my $aln (@sortedproj){
	    if($aln->[1] eq $seqname){
		my ($startcol,$stopcol) = AlignmentTree::coordstocolumn($atree->{_alignments}->{$aln->[0]}->[0],$seqname,$aln->[2],$aln->[3],1);
		my $sv = &reportVariants($atree,$db,$aln->[0],$seqname,$startcol,$stopcol,$nearestseq);
		foreach my $v (@$sv){
		    if($v->[4] != 0){
			print "#FSVAR $seqname ",join(',',@$v),"\n" if($debug);
			push @$fsvars,$v;
			$netfs += $v->[4];
			#if(abs($netfs) > $FS_THRESHOLD){
			    #Short circuit
			    #return undef;
			#}
		    }
		    else{
			die;
		    }
		}
	    }
	}		
	my @coords = sort {$a->[0] <=> $b->[0]} (@$fsvars);
	my $pos = [];
	my @runs;
	my $indelstr1;
	my $indelstr2;
	my $last;
	my $start;
	my $lasttype;
	my $end;
	for(my $i=0;$i<@coords;$i++){
	    if($i==0){
		$start=$coords[$i]->[0];
		$lasttype=$coords[$i]->[4];
		$pos= [];
	    }
	    elsif(abs($last+1 - $coords[$i]->[0]) > 1 || $coords[$i]->[4] != $lasttype){
		#print "Adding $start,$last,$indelstr1,$indelstr2,$lasttype ",scalar(@$pos),"\n";
		push @runs,[$start,$last,$indelstr1,$indelstr2,$lasttype,$pos];
		$indelstr1="";
		$indelstr2="";
		$start=$coords[$i]->[0];
		$pos=[];
	    }
	    $last=$coords[$i]->[0];
	    $lasttype=$coords[$i]->[4];
	    $indelstr1.=$coords[$i]->[1];
	    $indelstr2.=$coords[$i]->[2];
	    push @$pos,$coords[$i]->[0];
	}
	if($last){
	    #print "Adding_post $start,$last,$indelstr1,$indelstr2,$lasttype ",scalar(@$pos),"\n";
	    push @runs,[$start,$last,$indelstr1,$indelstr2,$lasttype,$pos];
	}
	my $ispmark=0;
	foreach my $r (@runs){
	    $ispmark = ($r->[2] =~ /$PMARK_SPACER/) ? 1 : 0;
	    $ispmark = ($r->[3] =~ /$PMARK_SPACER/) ? 1 : 0 if(!$ispmark);
	    last if($ispmark);
	}


	#Remove runs that are multiple of 3
	my @fsruns;
	foreach my $r (@runs){
	    die if(length($r->[2]) != length($r->[3]));
	    if($ispmark==1){
		push @fsruns,$r;
	    }else{#if(length($r->[2])%3!=0 || (length($r->[2]) < $FSLEN_THRESHOLD)){
		push @fsruns,$r;
	    }
	}
	
	if($verbose){
	    print "#FS num_runs ",scalar(@fsruns),"\n";
	    foreach my $r (@fsruns){
		print "[$r->[0]-$r->[1] $r->[2]:$r->[3]] $r->[4] ",scalar(@{$r->[5]}),"\n";
	    }
	    print "\n";
	}
	return (\@fsruns,$netfs);
    }
}

##############################
#Report ORFs on aligned sequences that are unannotated
#
#Method:

#For all aligned segments that do not contain any annotated ORF
#Attempt to use annotated and aligned start codons from other genomes in the cluster
#to call new ORFs
sub findnewORFs{
    my($db,$atree,$mappedorgs,$mappedgenes,$codons) = @_;
    #Consider all possible aligned starts in the cluster
    my $allcodons = {};
    foreach my $seq (keys %{$codons->{'starts'}}){
	foreach my $codon (keys %{$codons->{'starts'}->{$seq}}){
	    $allcodons->{$codon}++;
	}
    }
    my $noorfseqs = {};
    #Foreach codon, attempt to find ORFs if none annotated above cutoffs
    print "#Total number of possible codons ",scalar(keys %$allcodons),"\n" if($debug);;
    foreach my $codon (keys %$allcodons){
	my($col,$aln) = split(/\$CODON_DELIM_REGEX/,$codon);
	my $alignedseqs  = $atree->{_alignments}->{$aln}->[0]; #get seqs for $lan
	foreach my $alnseq (@$alignedseqs){
	    my $seq = $alnseq->[0];
	    #Check if sequence already has a mapped gene
	    if(! exists $mappedorgs->{$seq}){
		print "#No ORFs on seq $seq\n" if($debug);;
		$noorfseqs->{$seq}++;
	    }
	}
    }
    my $seq_attrs = {};
    print "#Looking for new starts in new orfs\n" if($debug);;
    &checkStarts($db,$codons,[keys %$noorfseqs],$seq_attrs,1);
    return $seq_attrs;
}

sub printExtAlts{
#Longest row in Green
#Frameshifts in Red
#Start fully consistent
#Alt,Num,Annotated,Length
    


}

sub printExtAlignments{
#Use CSS classes for each codon, highlight in color
#Use CSS classes for each gene

#getRowClass
#grid.getView().getRowClass = function(record, index){
#return (record.data.change<0.7 ? (record.data.change<0.5 ? (record.data.change<0.2 ? 'red-row' : 'green-row') : 'blue-row') : '');
#};

    var fDataTpl = new Ext.XTemplate(
        '<tpl for=".">',
            '<div>',
                '<pre class="x-fixed">{element}</pre>',
            '</div>',
        '</tpl>'
    );

#From http://www.sencha.com/blog/2010/07/13/a-side-by-side-diff-viewer-built-with-ext-js/
#        // Obtain reference to HTML templates
#        lineTpl = Ext.ux.CodeViewer.lineTpl,
#        emptyLineTpl = Ext.ux.CodeViewer.emptyLineTpl,

#        // Create a "pre" tag to hold the code
#        pre = this.el.createChild({tag: 'pre'}),

#        var el = lineTpl.append(pre, [i+1, this.highlightLine(lines[i])]);
# Ext.fly(el).addClass('ux-codeViewer-modified');


}




sub printExtJSCluster{

    my($cluster_id,$clusterref) = @_;
    #List cluster members
    my @clustergrid;

#Longest row in Green
#Frameshifts in Red
#Start fully consistent
#Alt,Num,Annotated,Length
    my @altgrid;
    foreach my $alt (keys %{$clusterref->{'alts'}}){
	my $isfcon = (exists $clusterref->{'alts'}->{$alt}->{'fcon'}) ? 1 : 0;
	my $ismax = (exists $clusterref->{'alts'}->{$alt}->{'maxlen'}) ? 1 : 0;
	push @altgrid,["'".$clusterref->{'alts'}->{$alt}->{'name'}."'",
		       $clusterref->{'alts'}->{$alt}->{'gfreq'},
		       $clusterref->{'alts'}->{$alt}->{'afreq'},
		       $clusterref->{'alts'}->{$alt}->{'len'},
		       scalar(keys %{$clusterref->{'alts'}->{$alt}->{'neworfs'}}),
		       scalar(keys %{$clusterref->{'alts'}->{$alt}->{'fs'}}),
		       $isfcon,
		       $ismax];
		       
    }
    foreach my $g (keys %{$clusterref->{'orgs'}}){
	my @codoninfo;
    	foreach my $alt (keys %{$clusterref->{'alts'}}){
	    if(exists $clusterref->{'codons'}->{$alt}->{'features'}){
		push @codoninfo,$clusterref->{'alts'}->{$alt}->{'name'};
	    }
	}
	my $gref = $clusterref->{'orgs'}->{$g};
	push @clustergrid,["'CLUSTER_".$cluster_id."'",
			   "'".$g."'",
			   "'".join(',',@{$gref->{'genes'}})."'",
			   "'".join(',',@{$gref->{'cov'}})."'",
			   "'".join(',',@{$gref->{'pid'}})."'",
			   $gref->{'fmin'},
			   $gref->{'fmax'},
			   $gref->{'len'},
			   "'".join(',',@{$gref->{'orient'}})."'",
			   "'".join(',',@codoninfo)."'",
			   "'".$gref->{'desc'}."'"
			   ];
    }
    #List edits

    #Show alignment

    if($htmlout){
	#Link to prev and next cluster    
	my $jsfh;
	my $htmlfh;
	my $htmlrelpath = basename("$options{'prefix'}cluster_${cluster_id}.html");
	my $jsrelpath = basename("$options{'prefix'}cluster_${cluster_id}.js");
	open $jsfh,"+>$options{'prefix'}cluster_${cluster_id}.js";
	open $htmlfh,"+>$options{'prefix'}cluster_${cluster_id}.html";
	
	print $htmlfh <<_CLUSTERHTMLHEADER;
	
	<html>
	    <head>
	    <title>Cluster $cluster_id</title>
	    <link rel="stylesheet" type="text/css" href="http://dev.sencha.com/deploy/dev/resources/css/ext-all.css" />
	    <script type="text/javascript" src="http://dev.sencha.com/deploy/dev/adapter/ext/ext-base.js"></script>
	    <script type="text/javascript" src="http://dev.sencha.com/deploy/dev/ext-all-debug.js"></script>
	    </head>
	    <body>
	    <script type="text/javascript" src="$jsrelpath"></script>
	    <div id="my-div" class="x-hidden">
            <pre>

_CLUSTERHTMLHEADER
;
	
	
	print $htmlfh `cat $options{'prefix'}cluster_${cluster_id}.aln.out`;
	
	print $htmlfh <<_CLUSTERHTMLFOOTER;

	    </pre>
	    
</body>
</html>

_CLUSTERHTMLFOOTER
;
	

    print $jsfh <<_CLUSTERJSHEADER;
    function renderGeneURL(val){
      return '<a href="javascript:document.getElementById(\\''+val+'\\').scrollIntoView(true);">'+val+'</a>';
    }
    Ext.onReady(function(){
	
	Ext.QuickTips.init();
	
	var xg = Ext.grid;
		
	var featstore = new Ext.data.ArrayStore({
	  fields: [
		   {name: 'cluster'},
		   {name: 'genome'},
		   {name: 'name'},
		   {name: 'coverage', type: 'float'},
		   {name: 'identity', type: 'float'},
		   {name: 'fmin', type: 'float'},
		   {name: 'fmax', type: 'float'},
		   {name: 'length', type: 'float'},
		   {name: 'strand'},
		   {name: 'codon_pairs'},
		   {name: 'desc'}
		   ]
	});
	
	var altstore = new Ext.data.ArrayStore({
	  fields: [
		   {name: 'name'},
		   {name: 'gfreq'},
		   {name: 'afreq'},
		   {name: 'len'},
		   {name: 'neworfs'},
		   {name: 'fs'},
		   {name: 'isfcon'},
		   {name: 'ismax'},
		   ]
	});
	
	featstore.loadData(xg.clusterData);
	altstore.loadData(xg.altData);

	var alntext = new Ext.Panel({
	    'id':'alntext',
	    'title':'Alignment detail',
	    'region':'south',
	  split:true,
	  height:300, 
	  collapsible: true,
	  autoScroll:true,
	  contentEl:'my-div'
	  });
	
	var featgrid = new xg.GridPanel({
	  store: featstore,
	  columns: [
		    {id:'cluster',header: "Cluster", width: 70, sortable: true, dataIndex: 'cluster'},
		    {header: "Feature", width: 100, sortable: true, dataIndex: 'name'},
		    {header: "Genome", width: 100, sortable: true, dataIndex: 'genome'},
		    {header: "fmin", width: 50, sortable: true, dataIndex: 'fmin'},
		    {header: "fmax", width: 50, sortable: true, dataIndex: 'fmax'},
		    {header: "strand", width: 30, sortable: true, dataIndex: 'strand'},
		    {header: "Len", width: 50, sortable: true, dataIndex: 'length'},
		    {header: "Coverage", width: 50, sortable: true, dataIndex: 'coverage'},
		    {header: "Identity", width: 50, sortable: true, dataIndex: 'identity'},
		    {header: "Alt ORFs", width: 200, sortable: true, dataIndex: 'codon_pairs'},
		    {header: "Description", width: 500, sortable: true, dataIndex: 'desc'},
		    ],	      
	    viewConfig: {
		    forceFit: true},
		    frame: true,
		    animCollapse: false,
		    title: 'Cluster ${cluster_id} annotation summary',
		    iconCls: 'icon-grid',
		    fbar  : ['->', {
		    text:'Save as text',
				     handler : null
			}],
		    columnWidth: .6,
		    flex:1
		});

	var editgrid = new xg.GridPanel({
		store: altstore,
		    columns: [
			      {id:'name',header: "ORF", width: 70, sortable: true, dataIndex: 'name'},
			      {header: "Aligned Freq", width: 50, sortable: true, dataIndex: 'gfreq'},
			      {header: "Annotated Freq", width: 50, sortable: true, dataIndex: 'afreq'},
			      {header: "Len", width: 50, sortable: true, dataIndex: 'len'},
			      {header: "# Missing", width: 50, sortable: true, dataIndex: 'neworfs'},
			      {header: "# FS", width: 50, sortable: true, dataIndex: 'fs'},
			      {header: "isfcon", width: 50, sortable: true, dataIndex: 'isfcon'},
			      {header: "ismax", width: 50, sortable: true, dataIndex: 'ismax'},
			      ],	      
		    split:true,
		    frame: true,
		    collapsible: true,
		    animCollapse: false,
		    title: 'Cluster 1 edit summary',
		    iconCls: 'icon-grid',
		    columnWidth: .4,
		    flex:1
		    });
	var viewport = new Ext.Viewport({
	  layout:'border',
	    //defaults: {autoScroll:true,height:500},
	  items: [
		  new Ext.Panel({
		    layout:'fit',
		    region:'center',
		    items: [
			    new Ext.Panel({
			      layout:'hbox',
			      layoutConfig: {
				  align : 'stretch',
				  pack  : 'start',
			      },
							    region:'center',
			      items: [ featgrid,editgrid]
			      })
						]
					    }),
		  alntext
		  ]
	      });
	viewport.doLayout();

    });

_CLUSTERJSHEADER
    ;
 
    print $jsfh "Ext.grid.clusterData = [";
    foreach my $c (@clustergrid){
	print $jsfh "[",join(',',@$c),"],\n";
    }
    print $jsfh "];\n";

    print $jsfh "Ext.grid.altData = [";
    foreach my $c (@altgrid){
	print $jsfh "[",join(',',@$c),"],\n";
    }
    print $jsfh "];\n";

    close $jsfh;
    close $htmlfh;
    }

}

sub printExtJS{
    my($clusters) = @_;
    my @clustergrid;
    
    #Summary is cluster_id,#genomes,#genes,class,

    foreach my $cluster_id (keys %$clusters){
	push @clustergrid,[$cluster_id,$clusters->{$cluster_id}->{'num_feats'},$clusters->{$cluster_id}->{'num_genomes'},"'".$clusters->{$cluster_id}->{'classes'}."'"];
	&printExtJSCluster($cluster_id,$clusters->{$cluster_id});
    }

    my $jsfh;
    my $htmlfh;
    my $jsrelpath = basename("$options{'prefix'}main.js");
    my $relpath = basename("$options{'prefix'}");
    open $jsfh,"+>$options{'prefix'}main.js";
    open $htmlfh,"+>$options{'prefix'}index.html";

    print $htmlfh <<_HTMLHEADER;

<html>
<head>
<title>Mugsy-Annotator Report</title>
<link rel="stylesheet" type="text/css" href="http://dev.sencha.com/deploy/dev/resources/css/ext-all.css" />
<script type="text/javascript" src="http://dev.sencha.com/deploy/dev/adapter/ext/ext-base.js"></script>
<script type="text/javascript" src="http://dev.sencha.com/deploy/dev/ext-all-debug.js"></script>
</head>
<body>
<script type="text/javascript" src="$jsrelpath"></script>

</body>
</html>

_HTMLHEADER
;
    

    print $jsfh <<_MAINJSHEADER;
	

    function renderClusterURL(val){
	return '<a href="${relpath}cluster_'+val+'.html">CLUSTER_'+val+'</a>';
    }

    Ext.onReady(function(){

	Ext.QuickTips.init();
	
	var xg = Ext.grid;
	
	
	// shared reader
	    var reader = new Ext.data.ArrayReader({}, [
						       {name: 'cluster_id'},
						       {name: 'num_feats'},
						       {name: 'num_genomes'},
						       {name: 'quality_class'}
						       ]);
	var store = new Ext.data.GroupingStore({
	  reader: reader,
	  data: xg.summaryData,
	  sortInfo:{field: 'cluster_id', direction: "ASC"},
	  groupField:'quality_class'
        });
	
	var grid = new xg.GridPanel({
	  store: store,
	  columns: [
		    {id:'Cluster',header: "Cluster", width: 10, sortable: true, dataIndex: 'cluster_id', renderer:renderClusterURL},
		    {header: "Features", width: 10, sortable: true, dataIndex: 'num_feats'},
		    {header: "Genomes", width: 10, sortable: true, dataIndex: 'num_genomes'},
		    {header: "Class", width: 20, sortable: true, dataIndex: 'quality_class'},
		    ],
	      
        view: new Ext.grid.GroupingView({
            forceFit:true,
            groupTextTpl: '{text} ({[values.rs.length]} {[values.rs.length > 1 ? "Items" : "Item"]})'
        }),

        frame:true,
        width: 700,
        height: 450,
        collapsible: true,
        animCollapse: false,
        title: 'Annotation summary',
        iconCls: 'icon-grid',
        fbar  : ['->', {
            text:'Clear Grouping',
            iconCls: 'icon-clear-group',
            handler : function(){
                store.clearGrouping();
            }
        }],
        renderTo: document.body
    });
});

_MAINJSHEADER
    ;


    print $jsfh "Ext.grid.summaryData = [";

    foreach my $c (@clustergrid){
	print $jsfh "[",join(',',@$c),"],\n";
    }
    print $jsfh "];\n";

    close $jsfh;
    close $htmlfh;

}

