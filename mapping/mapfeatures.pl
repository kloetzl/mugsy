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
-checkFrameShifts doesn't work on - strand

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
#use lib '/usr/local/projects/angiuoli/developer/sangiuoli/mugsy/trunk/mapping';

use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

#Bioperl is used only for translation machinery
use Bio::Perl;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::Tools::CodonTable;
use Bio::Seq::EncodedSeq;
#use Bio::LiveSeq::Mutation; tried this but couldn't get to work properly
#Default cutoffs

use AlignmentTree;

my %options;
my $results = GetOptions (\%options, 
			  'prefix=s',
			  'map_file=s',
			  'coverage|c=s',
			  'query_coverage|q=s',
			  'identity|i=s',
			  'sortkeys=s', #
			  'reportedits=s', #number of edits to report
			  'maxchange=s', #max allowable %length changes
			  'prefix=s', #Generate output reports with file prefix
			  'cogformat=s',
			  'printalignments',
			  'skipaltstarts',
			  'skipneworfs',
			  'showframeshifts',
			  #Gene calling options
			  'minorflen=s',
			  'maxorflen=s',
			  'verbose|v',
			  'debug|d=s') || pod2usage(-verbose => 1);

pod2usage(-verbose=>1) if($options{'help'});


my $coverage_cutoff = (exists $options{'coverage'}) ?  $options{'coverage'} : 0.5;
my $query_coverage_cutoff = (exists $options{'query_coverage'}) ?  $options{'query_coverage'} : 0;
my $pid_cutoff= (exists $options{'identity'}) ?  $options{'identity'} : 0.1;
print STDERR "Using coverage cutoff:$coverage_cutoff identity:$pid_cutoff query_coverage:$query_coverage_cutoff\n";

my $MAXORFLEN = (exists $options{'maxorflen'}) ? $options{'maxorflen'} : 30000; #in bp
my $MINORF= $options{'minorflen'} || 50; #in aa residues
my $FS_WINDOW = 200; #bp window upstream from pre-mature stop to check for frameshifts
my $FS_THRESHOLD = 3;
my $ORFLEN_MAXDELTA = 0.5; #do not consider possible codons that are less than X the length of the maximum annotated ORF

#Used for detecting contig boundaries
my $PMARK_SPACER = "NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN";

#Flag for checking consistent start,stop
#Assumes input features are genes
my $doconsistencychecks=0;

#Report new ORFs using aligned start codons
my $dofindneworfs = 1;
my $autocorrect=0;

#my $autofixunmapped=0; #TODO, use consistency checks to fix annotations of unmapped genes

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
my $printskipped=1;
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
my $verbose=$options{'verbose'};

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
my $codons = {};
my $classes_sum = {};
my $newclasses_sum = {};

#AlignmentTree is a interval tree that contains alignments between sequences
#and features on those sequences
my $atree = AlignmentTree::deserialize($ARGV[0]);
$atree->{_debug}=$debug;
#Read a white space delimited list of features to map from stdin

my %featlookup;
my $filetype;
my $fh;

if($options{'map_file'}) {
	open($fh, "<$options{'map_file'}") or die "Error in opening the file, $options{'map_file'}, $!\n";
} else {
	$fh = \*STDIN;
}

while(my $line=<$fh>){
    my($name,$seq,$fmin,$fmax,$orient,$polyid,$geneid,$annotations);
    chomp $line;
    if($line =~ /\#gff-version 3/){
	$filetype = 'gff3';
	$featlookup{'gene'}++;
	$featlookup{'pseudogene'}++;
    }
    if($line !~ /^\#/){
	if($filetype eq 'gff3'){
	    #GFF
	    my @elts = split(/\t/,$line);
	    if(scalar(@elts)==9){
		if(exists $featlookup{lc($elts[2])}){
		    my %attrs = map {split(/=/)} split(/;/,$elts[8]);
		    if(exists $attrs{'locus_tag'}){
			$name = $attrs{'locus_tag'};
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
		    my $cdsname;
		    if(exists $attrs{'locus_tag'}){
			$cdsname = $attrs{'locus_tag'};
		    }
		    if(exists $attrs{'product'}){
			$annotations .= $attrs{'product'};
			
		    }
		    if(exists $features->{$cdsname}){
			$features->{$cdsname}->[11] = $annotations;
		    }
		}
	    }
	    else{
		print "Skipping $line\n";
	    }
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
	    die "Bad orient $orient. $line" if($orient ne '+' && $orient ne '-');
	    $features->{$name} = [$seq,$fmin,$fmax,$fmax-$fmin,$orient,$polyid,$geneid];
	    #[7]-[10] reserved for start,stop codon info
	    $features->{$name}->[11] = $annotations;
	}
    }
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
my $neworfcount = 0;
my $adjustedorfs = 0;

#List of newly called ORFs
my $neworfs = {};
#and the annotated ORFs they replace
my $subsumed = {};

#Map of feature => organism
my $feat2organism = {};

my $db;

if(! $options{'prefix'}){
    $options{'prefix'} = "mugsyant.$$";
}

open CFILE, "+>$options{'prefix'}clusters.out";

if(-f "$ARGV[1]"){
    print STDERR "Using FASTA file $ARGV[1]\n";
    $db = Bio::DB::Fasta->new($ARGV[1],'-reindex'=>1); 
    my @ids = $db->ids();
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

    #As the algorithm progresses, features are mapped and removed from consideration
    #Consider genes that remain unmapped or 
    #remain covered by <= cutoff% of length in alignments already considered
    if(!exists $mapped->{$query} && !exists $deleted->{$query}){

	#Start a new cluster based on the query gene.  Set a new
	#cluster id; each cluster can also be identified by the query
	#gene ($query)
	$cluster_id++;
	
	my($mappedorgs,$mappedgenes,$unmappedorgs,$unmappedgenes) = &buildCluster($atree,$query);
	print "MAPPED Num_orgs:",scalar(keys %$mappedorgs)," Num_genes:",scalar(keys %$mappedgenes)," UNMAPPED Num_orgs:",scalar(keys %$unmappedorgs)," Num_genes:",scalar(keys %$unmappedgenes),"\n" if($debug);
	die "Less than 2 mapped sequences" if(scalar(keys %$mappedgenes)>1 && scalar(keys %$mappedorgs)<=0);
	die "No mapped genes" if(scalar(keys%$mappedgenes)<1);

	#Mark inconsistencies in the cluster and save codons
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
	if($dofindneworfs && !$options{'skipneworfs'}){
	    $new_orfs = &findnewORFs($db,$atree,$mappedorgs,$mappedgenes,$codons);	    
	}
	

	if((scalar(keys %$mappedgenes)>1 && scalar(keys %$mappedorgs)>1)){
	    print "#Cluster WGA$cluster_id\n" if($debug);;
	    if($doconsistencychecks){
		#Determine annotation inconsistencies
		#(1) Check for genes called on the wrong strand
		#foreach my $feat_name (keys %$mappedgenes){
		#if($mappedgenes->{$feat_name}->{'morient'}!=0){
		#print "#Class O1. Mis-matched orientation on $feat_name\n" if($debug);;
		#$feat_attrs->{$feat_name}->{'classO1'} = $mappedgenes->{$feat_name}->{'morient'};
	        #}
	        #}
		
		#A consistent cluster has both consistent start and stop codons (labeled CS1,CE1)
		my $consistent=(exists $cluster_attrs->{'CS1'} && exists $cluster_attrs->{'CE1'});
		
		if(!$consistent){
		    #(2) Attempt to adjust start sites
		    #Using aligned codons as possible alternative start sites
		    #Print out alternative (fixed) gene coordinates		    
		    print "#Looking for new starts\n" if($debug);;
		    if(!defined $options{'skipaltstarts'}){
			my $altorfs = &checkStarts($db,$codons,[keys %$mappedorgs],$seq_attrs);
			#Report alternative start codons in decreasing order of use
			foreach my $aorf (@$altorfs){
			    my($seqname,$fmin,$fmax,$orforient,$orfstartcodon,$orfstopcodon) = @_;
			    my $frame;
			    my $dist;
			    if($orforient eq '-'){
				$frame=(($db->get_Seq_by_id($seqname)->length()-$fmax)%3)*-1;
				$dist=0;#$end-$althash->{'end'};
			    }
			    else{
				$frame=$fmin%3;
				$dist=0;#$start-$althash->{'start'};
			    }
			    
			    if(!exists $seq_attrs->{$seqname}){
				$seq_attrs->{$seqname} = [];
			    }
			    if($freq_altstops){ 
				#check if stop is consistent with at least one other stop
				#if($althash->{'stopfreq'}>1){
				    #push @{$seq_attrs->{$seqname}},"alt_start=$althash->{'start'}-$althash->{'end'},orient:$althash->{'orient'},len:".($althash->{'end'}-$althash->{'start'}).",pairfreq:$codons->{'pairs'}->{$althash->{'startcodon'}.':'.$althash->{'stopcodon'}}->{'gfreq'},startfreq:$althash->{'startfreq'},stopfreq:$althash->{'stopfreq'},frame:$frame,dist:$dist,startcodon:$althash->{'startcodon'},stopcodon:$althash->{'stopcodon'};";
			    #}
			    }
			    else{
				#push @{$seq_attrs->{$seqname}},"alt_start=$althash->{'start'}-$althash->{'end'},orient:$althash->{'orient'},len:".($althash->{'end'}-$althash->{'start'}).",pairfreq:$codons->{'pairs'}->{$althash->{'startcodon'}.':'.$althash->{'stopcodon'}}->{'gfreq'},startfreq:$starts->{$althash->{'startcodon'}},stopfreq:$stops->{$althash->{'stopcodon'}},frame:$frame,dist:$dist,startcodon:$althash->{'startcodon'},stopcodon:$althash->{'stopcodon'};";
			    }
			}
		    }
		    #foreach my $feat_name (keys %$mappedgenes){
			#$feat_attrs->{$feat_name}->{'pairfreq='.$codons->{'pairs'}->{"$codons->{'featstarts'}->{$feat_name}"."$codons->{'featstops'}->{$feat_name}"}}++;
			#$feat_attrs->{$feat_name}->{'pairfreq'}=$codons->{'pairs'}->{"$codons->{'featstarts'}->{$feat_name}".':'."$codons->{'featstops'}->{$feat_name}"}->{'gfreq'};
		    #}
		    #(3) Attempt to resolve inconsistencies by frameshifting
		    #Using aligned codons as start sites and indels adjacent to stops
		    #as possible locations of the frameshift
		    #Long ORFs may indicate a sequencing error or authentic frameshift
		    #Reports 
		    #-frameshift location and indel
		    #-new ORF,translation
		    #-any deleted genes
		    if(defined $options{'showframeshifts'}){
			print "#Looking for frameshifts\n" if($debug);;
			#&checkFrameshifts($db,$mappedorgs,$mappedgenes,$codons,$seq_attrs);
			&checkFrameshifts_new($db,$codons,$mappedorgs,$seq_attrs);
		    }
		}
	    }

	    
	    
	
	    #We have a good cluster, save it
	    #Save the cov,pid in master list of mapped genes

	    my $totallen=0;
	    my $maxlen=0;
	    foreach my $feat_name (keys %$mappedgenes){
		die "Feature $feat_name already mapped" if(exists $mapped->{$feat_name});
		$mapped->{$feat_name}->{'cov'}=$mappedgenes->{$feat_name}->{'cov'}/$features->{$feat_name}->[3];
		$mapped->{$feat_name}->{'pid'}=$mappedgenes->{$feat_name}->{'pid'}/$mappedgenes->{$feat_name}->{'len'};
		$totallen += $features->{$feat_name}->[3];
		$maxlen = ($features->{$feat_name}->[3] > $maxlen) ? $features->{$feat_name}->[3] : $maxlen;
		delete $unmapped->{$feat_name};
	    }
	    my $avglen=$totallen/(scalar keys %$mappedgenes);
	    my $fshifts=0;
	    my $classesstr;
	    $debug=1;
	    if(!defined $options{'reportedits'} || $options{'reportedits'} > 0){
		#Save aligned and annotated codon frequency
		foreach my $p (keys %{$codons->{'pairs'}}){
		    print "#Analyzing codon pair $p\n" if($debug);
		    my($startcodon,$stopcodon) = split(/:/,$p);
		    foreach my $seqname (keys %$mappedorgs,keys %$unmappedorgs){
			print "#Sequence $seqname\n" if($debug);
			#if this is the annotated pair
			if(exists $codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} && $codons->{'pairs'}->{$p}->{'orgs'}->{$seqname}->[3]==1){
			    #Do nothing, already annotated
			    $codons->{'pairs'}->{$p}->{'features'}->{$seqname} = $mappedorgs->{$seqname}->{'features'};	
			    print "#annotated\n" if($debug);
			}
			else{
			    #check if this is an ORF in $seqname
			    my($fmin,$fmax,$orient) = &findCoords($atree,$seqname,$startcodon,$stopcodon);
			    if(defined $fmin && defined $fmax && defined $orient && $fmax>$fmin){
				die "$atree,$seqname,$startcodon,$stopcodon" if(! defined $fmin || ! defined $fmax);
				if(&isORF($db,$seqname,$fmin,$fmax,$orient)){
				    if(exists $unmappedorgs->{$seqname}){
					#Requires a new ORF in an unannotated region
					$codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,-1];
					print "#neworf $p $seqname ",$fmax-$fmin,"\n" if($debug);

				    }
				    else{
					#Requires a new ORF that is different than currently annotated
					$codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,0];
					if(exists $mappedorgs->{$seqname}){
					    $codons->{'pairs'}->{$p}->{'features'}->{$seqname} = $mappedorgs->{$seqname}->{'features'};	
					    print "#altorf\n" if($debug);
					}
					else{
					    print "#altorf, prev did not pass cutoffs\n" if($debug);
					}
				    }
				}
				else{
				    #Look for possible frameshifts if there are either
				    #a) Multiple annotated ORFs in this region
				    #b) An upstream stop codon
				    my $feat_name;
				    my $annotatedstop;
				    if(exists $mappedorgs->{$seqname} && scalar(keys %{$mappedorgs->{$seqname}->{'features'}}) == 1){
					$feat_name = [keys %{$mappedorgs->{$seqname}->{'features'}}]->[0];
					$annotatedstop = $features->{$feat_name}->[7] . $CODON_DELIM . $features->{$feat_name}->[8];
				    }
				    if(exists $unmappedorgs->{$seqname}
				       || 
				       (exists $mappedorgs->{$seqname}
					&&
					(scalar(keys %{$mappedorgs->{$seqname}->{'features'}}) > 1
					 ||
					 $stopcodon ne $annotatedstop
					 )
					)){
					die "$seqname found in both mapped and unmapped org lists" if(exists $unmappedorgs->{$seqname} && exists $mappedorgs->{$seqname});
					print "#Considering FS $stopcodon ne $annotatedstop for $feat_name on $seqname\n" if($debug && exists $mappedorgs->{$seqname});
					#Find most similar sequence that has this ORF
					my @neighborseqs = keys %{$codons->{'pairs'}->{$p}->{'orgs'}};
					my($nearestseq,$indels,$fs) = &findNearestNeighbor($atree,$seqname,\@neighborseqs,$startcodon,$stopcodon);
					print "#Using $nearestseq as nearest neighbor to $seqname\n" if($debug);
					#Look for frameshifting mutations in $seqname
					my($fs,$netfs) = &reportFrameShifts($atree,$db,$seqname,$nearestseq,$startcodon,$stopcodon);
					if(ref $fs){
					    print "#Possible ORF with frameshift #indels:",scalar(@$fs)," net:$netfs\n" if($debug);
						if(abs($netfs) <= $FS_THRESHOLD){
						    print "#Adding frameshift net:",scalar(@$fs)," $netfs\n" if($debug);
						    $codons->{'pairs'}->{$p}->{'orgs'}->{$seqname} = [$fmin,$fmax,$orient,0,$fs];
						    if(exists $mappedorgs->{$seqname}){
							$codons->{'pairs'}->{$p}->{'features'}->{$seqname} = $mappedorgs->{$seqname}->{'features'};
						    }
						}
					}
				    }
				}
			    }
			}
		    }
		}
		foreach my $p (keys %{$codons->{'pairs'}}){
		    foreach my $org (keys %{$codons->{'pairs'}->{$p}->{'orgs'}}){
			print "#CODONPAIR ",join(',',@{$codons->{'pairs'}->{$p}->{'orgs'}->{$org}}),"\n" if($verbose);
			$codons->{'pairs'}->{$p}->{'gfreq'}++;
			$codons->{'pairs'}->{$p}->{'afreq'}++ if($codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[3]==1); #inc only if annotated
			$codons->{'pairs'}->{$p}->{'length'}+=($codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[1] - $codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[0]);
			if(ref $codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[4]){
			    $codons->{'pairs'}->{$p}->{'fsvars'}+=1;
			    print "#FS ",join(',',@{$codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[4]}),"\n" if($verbose);
			}
			$codons->{'pairs'}->{$p}->{'neworfs'}+=1 if($codons->{'pairs'}->{$p}->{'orgs'}->{$org}->[3]<0);
		    }
		$codons->{'pairs'}->{$p}->{'len'} = $codons->{'pairs'}->{$p}->{'length'}/$codons->{'pairs'}->{$p}->{'gfreq'} if($codons->{'pairs'}->{$p}->{'gfreq'} > 0);
		}

		$classesstr = join(';',sort {$a cmp $b} keys %{$cluster_attrs});
		#Suggest edits for inconsistently annotated clusters
		if($classesstr ne 'CE1;CS1'){
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
				
				my $deltafracmax = abs($codonlength-$maxlen)/$maxlen;
				if($deltafracmax < $ORFLEN_MAXDELTA){
				    print EFILE ">CLUSTER_$cluster_id $bestcodon\n";
				    my $newmappedorgs; 
				    my $newmappedgenes;
				    foreach my $org (keys %{$codons->{'pairs'}->{$bestcodon}->{'orgs'}}){
					if(scalar(keys %{$codons->{'pairs'}->{$bestcodon}->{'features'}->{$org}})>0){
					    foreach my $feat_name (keys %{$codons->{'pairs'}->{$bestcodon}->{'features'}->{$org}}){
						my $pred_feat = $codons->{'pairs'}->{$bestcodon}->{'orgs'}->{$org};
						if($pred_feat->[0] ne $features->{$feat_name}->[1] || 
						   $pred_feat->[1] ne $features->{$feat_name}->[2]){
						    my $fs = (defined $codons->{'pairs'}->{$bestcodon}->{'orgs'}->{$org}->[4]) ? "F" : "";
						    print EFILE "$feat_name\t$org\t$pred_feat->[0]\t$pred_feat->[1]\t",($pred_feat->[1] - $pred_feat->[0]),"\t$pred_feat->[2]\t$fs\n";
						}
						$newmappedorgs->{$org}->{'features'}->{$feat_name}++;
						$newmappedgenes->{$feat_name}->{'fmin'} = $pred_feat->[0];
						$newmappedgenes->{$feat_name}->{'fmax'} = $pred_feat->[1];
						$newmappedgenes->{$feat_name}->{'len'} = $pred_feat->[1] - $pred_feat->[0];
						$newmappedgenes->{$feat_name}->{'relorient'} = $pred_feat->[2];
					    }
					}
					else{
					    $cluster_attrs->{'CN1'}++;
					    my $feat_name = "NEWORF_CLUSTER$cluster_id";
					    my $pred_feat = $codons->{'pairs'}->{$bestcodon}->{'orgs'}->{$org};
					    print EFILE "$feat_name\t$org\t$pred_feat->[0]\t$pred_feat->[1]\t",($pred_feat->[1] - $pred_feat->[0]),"\t$pred_feat->[2]\n";
					    $newmappedorgs->{$org}->{'features'}->{$feat_name}++;
					    $newmappedgenes->{$feat_name}->{'fmin'} = $pred_feat->[0];
					    $newmappedgenes->{$feat_name}->{'fmax'} = $pred_feat->[1];
					    $newmappedgenes->{$feat_name}->{'len'} = $pred_feat->[1] - $pred_feat->[0];
					    $newmappedgenes->{$feat_name}->{'relorient'} = $pred_feat->[2];
					}
				    }
				    my $newclassesstr = $codons->{'pairs'}->{$bestcodon}->{'cluster_attrs'};
				    print "#BEST CODON $bestcodon $codons->{'pairs'}->{$bestcodon}->{'gfreq'} max_annotated_len:$maxlen delta_len_max:$deltafracmax\n"; 
				    print "#EDITTBL CLUSTER_$cluster_id $bestcodon\t$codons->{'pairs'}->{$bestcodon}->{'gfreq'}";
				    if($codons->{'pairs'}->{$bestcodon}->{'fsvars'} > 0){
					print "(F:$codons->{'pairs'}->{$bestcodon}->{'fsvars'})";
				    }
				    if($codons->{'pairs'}->{$bestcodon}->{'neworfs'} > 0){
					print "(N:$codons->{'pairs'}->{$bestcodon}->{'neworfs'})";
				    }
				    
				    print "\t$codons->{'pairs'}->{$bestcodon}->{'afreq'}\t$codons->{'pairs'}->{$bestcodon}->{'len'}\t$codons->{'pairs'}->{$bestcodon}->{'fshifts'}\t$codons->{'pairs'}->{$bestcodon}->{'neworfs'}\t";
				    if($codonlength eq $maxlen){
					print "#MAXLENEDIT ";
				    }
				    if($codons->{'pairs'}->{$bestcodon}->{'gfreq'} eq scalar(keys %$mappedorgs)){
					print "#FCONSISTENT ";
				    }
				    print "#$classesstr orig\n";
				    #print "#$newclassesstr new\n";
				    #&reportCluster($query,$newmappedorgs,$newmappedgenes,{},$feat_attrs,$cluster_attrs,$seq_attrs);
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
	    $classesstr = join(';',sort {$a cmp $b} keys %{$cluster_attrs});
	    #Print cluster
	    $debug=0;
	    &reportCluster($query,$mappedorgs,$mappedgenes,$unmappedgenes,$feat_attrs,$cluster_attrs,$seq_attrs,$new_orfs);

	    $classes_sum->{$classesstr}->{'ngenes'} +=scalar(keys %$mappedgenes);
	    $classes_sum->{$classesstr}->{'nclusters'}++;

	    $validcluster++;
	    if($COGoutputformat){}
	    else{
		print "#VALID\tCLUSTER_$cluster_id\tNum_organisms=",scalar(keys %$mappedorgs)+1,
		"\tNum_genes=",scalar(keys %$mappedgenes),"\n" if($debug);;
	    }
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
	    if($COGoutputformat){}
	    else{
		print "#SKIPPED\t$query\tWGA$cluster_id\tNum_organisms=",scalar(keys %$mappedorgs),
		"\tNum_genes=",scalar(keys %$unmappedgenes),"\n" if($debug);
	    }
            #Cluster does not pass cutoffs
	    #This cluster was skipped because it does not pass coverage cutoffs
	    #Optionally print
	    if($printskipped){
		if($COGoutputformat){
		}
		else{
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
	}
	

	foreach my $organism (keys %$new_orfs){
	    my $orfidx=0;
	    foreach my $alt (@{$new_orfs->{$organism}}){
		$neworfcount++;
	    }
	}

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



print "#NUM CLUSTERS $validcluster\n";

#Mark the remaining features as singletons categorized as
#1)not found in any alignments !exists mapped && !exists unmapped
#2)aligned but below cutoffs   !exists mapped && exists unmapped
$nomatches = &findSingletons($atree,$mapped,$unmapped,$subsumed);

#Calculate summary stats 
my $avgcov=0;
my $avgid=0;
my $mappedgenescount=0;

my $unmappedgenescount=0;
my $avgunmappedcov=0;
my $avgunmappedid=0;

my $nohit=0;

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

die if($nohit != scalar(keys %$nomatches));
   
#Print summary stats
print "\n\n\n";
print "Class legend\n";
print "C{S,E}1 - consistent start,stop\n";
print "C{S,E}2 - inconsistent start,stop\n";
print "C{S,E}3 - unaligned start,stop\n";
print "C{S,E}0 - start,stop is missing, runs off a contig end\n";
print "CM1 - multiple fragments spanned\n";
print "CO1 - mismatch orientation\n";
print "CN1 - aligned ORFs on unannotated genomes\n";

foreach my $cstr (sort {$classes_sum->{$b}->{'ngenes'} <=> $classes_sum->{$a}->{'ngenes'}} (keys %$classes_sum)){
    print "$cstr: num_genes:$classes_sum->{$cstr}->{'ngenes'} num_clusters:$classes_sum->{$cstr}->{'nclusters'}\n";
}
print "Number of clusters containing aligned features\n";
print "CLUSTERS: $validcluster\n";
print "Number aligned features mapped into clusters\n";
print "MAPPED: $mappedgenescount AVGCOV:",$avgcov/$mappedgenescount," AVGID:",$avgid/$mappedgenescount,"\n" if($mappedgenescount);
print "Number features with an overlapping alignment but are not mapped into clusters\n";
print "UNMAPPED: $unmappedgenescount AVGCOV:",$avgunmappedcov/$unmappedgenescount," AVGID:",$avgunmappedid/$mappedgenescount,"\n" if($unmappedgenescount && $mappedgenescount);
print "Number of features with no overlapping alignment\n";
print "NOHIT:$nohit\n";
print "Number of new ORFs that match an annotated start and stop codon\n";
print "NEWORF:$neworfcount\n";
print "POST CORRECTION alt starts, long orfs. Deleted ORFs ",scalar(keys %$deleted)," Adjusted ORFs:",$adjustedorfs,"\n";
foreach my $cstr (sort {$newclasses_sum->{$b}->{'ngenes'} <=> $newclasses_sum->{$a}->{'ngenes'}} (keys %$newclasses_sum)){
    print "$cstr: num_genes:$newclasses_sum->{$cstr}->{'ngenes'} num_clusters:$newclasses_sum->{$cstr}->{'nclusters'}\n";
}
print "POST CORRECTION alt starts, long orfs w/ frameshifts\n";

close CFILE;
close EFILE;
close DFILE;

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
    my $qcurrorg = '?';
    my $qcov = 0;
    my $qpid = 0;
    my $qalnfmin = undef;
    my $qalnfmax = undef;
    my $qfmin = $features->{$query}->[1];
    my $qfmax = $features->{$query}->[2];
    my $qfeatlen = $qfmax-$qfmin;
    my $qrelorient = 0;

    my $orient = $features->{$query}->[4];
    
    print "#MAPFEATURE Mapping $query $qseqname:$qfmin-$qfmax len:",$qfmax-$qfmin,"\n" if($debug);;
    
    #AlignmentTree::map() 
    #returns [0=alignment_name,1=seqid,2=align_start,3=align_stop,4=align_cov,5=feature_name,6=seqid,7=feature_cov,8=feature_pid]

    my @isect = $atree->map($qseqname,$qfmin,$qfmax,"alignment");
    
    #List of alignments that comprise the current cluster
    my $goodalignments = {};
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
    
    #First screen all overlapping alignments to ensure that they
    #include the query gene
    foreach my $r (@isect){
	my $feat_name = $r->[0];
	my $seqname = $r->[1];
	my $align_name = $r->[5];
	#Only consider WGA alignments (alignment name in $align_name) that span query (gene name in $feat_name)
	if($feat_name eq 'gene:'.$query){
	    print "#Mapped $feat_name $query $align_name\n" if($debug);
	    $goodalignments->{$align_name}++;
	    my $alignedseqs  = $atree->{_alignments}->{$align_name}->[0];
	    foreach my $seq (@$alignedseqs){
		die if(ref $seq->[0]);
		$allseqs->{$seq->[0]}++;
	    }
	}
    }
    
    if($COGoutputformat){}
    else{
	print "#QUERY=$query coords=$qfmin-$qfmax len=$features->{$query}->[3] strand=$features->{$query}->[4] Num_alignments=",scalar(keys %$goodalignments),"\n";
    }
    
    #Transform feat_name
    my @nisect;
    foreach my $r (@isect) {
	$r->[0] =~ s/gene\://;
	push @nisect,$r if(exists $features->{$r->[0]});
    }

    foreach my $r (
		   sort { $features->{$b->[0]}->[3] <=> $features->{$a->[0]}->[3] } #sort on feature length
		   
                          #sort {   #Sort alphanumeric on alignment_name, secondary on align_start
			  #    if($a->[0] eq $b->[0]){
			  #	   $a->[2] <=> $b->[2];
			  #    }
			  #    else{
			  #	   $a->[0] cmp $b->[0];
			  #    }
			  #} 
		   @nisect){
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
    $mappedgenes->{$query}->{'pid'} = $qpid;#TODO, the pid is wrong is wrong
    $mappedgenes->{$query}->{'len'} = $qfeatlen;
    $mappedgenes->{$query}->{'relorient'} = $qrelorient;
    $mappedgenes->{$query}->{'alignments'} = [keys %$goodalignments];
    $mappedorgs->{$qcurrorg}->{'features'}->{$query}++;
    $mappedorgs->{$qcurrorg}->{'qcov'} = $qalnfmax-$qalnfmin;


    #Set query coverage
    foreach my $feat_name (keys %$alnfeats){
	my $fmin = $features->{$feat_name}->[1];
	my $fmax = $features->{$feat_name}->[2];
	my $featlen = $fmax-$fmin;
	if($alnfeats->{$feat_name}->{'cov'}/$featlen >= $coverage_cutoff  && #%coverage over matching gene length
	   $alnfeats->{$feat_name}->{'pid'}/$alnfeats->{$feat_name}->{'len'} >= $pid_cutoff){ #%id over aligned length only
	    print "Summing query coverage feat_name $feat_name $feat2organism->{$feat_name} = $alnfeats->{$feat_name}->{'qcov'}. Current total $alnorgs->{$feat2organism->{$feat_name}}->{'qcov'}\n" if($debug);
	    $alnorgs->{$feat2organism->{$feat_name}}->{'qcov'} += $alnfeats->{$feat_name}->{'qcov'};
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

	if(($query_coverage_cutoff==0 || ($alnorgs->{$feat2organism->{$feat_name}}->{'qcov'}/$qfeatlen >= $query_coverage_cutoff)) &&
	   $alnfeats->{$feat_name}->{'cov'}/$featlen >= $coverage_cutoff  && #%coverage over matching gene length
	   $alnfeats->{$feat_name}->{'pid'}/$alnfeats->{$feat_name}->{'len'} >= $pid_cutoff){ #%id over aligned length onl
	    
	    print "PASSED\n" if($debug);
	    
	    #Check matching len is <= length of gene
	    $mappedorgs->{$feat2organism->{$feat_name}}->{'features'}->{$feat_name}++;
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
	}
	else{
	    print "BELOW\n" if($debug);
	    #Does not pass cutoffs	    
	    $unmappedgenes->{$feat_name}->{'cov'} = $alnfeats->{$feat_name}->{'cov'};
	    $unmappedgenes->{$feat_name}->{'fmin'} = $alnfeats->{$feat_name}->{'fmin'};
	    $unmappedgenes->{$feat_name}->{'fmax'} = $alnfeats->{$feat_name}->{'fmax'};
	    $unmappedgenes->{$feat_name}->{'pid'} = $alnfeats->{$feat_name}->{'pid'};
	    $unmappedgenes->{$feat_name}->{'len'} = $alnfeats->{$feat_name}->{'len'};
	    $unmappedgenes->{$feat_name}->{'relorient'} = $alnfeats->{$feat_name}->{'relorient'};

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
	    $unmappedorgs->{$feat2organism->{$feat_name}}->{'features'}->{$feat_name}++;
	    $unmappedorgs->{$feat2organism->{$feat_name}}->{'qcov'} = $alnorgs->{$feat2organism->{$feat_name}}->{'qcov'};
	}
    }
    return($mappedorgs,$mappedgenes,$unmappedorgs,$unmappedgenes);
}
    


#Classify consistency of annotations within a cluster
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
	my($startcodon,$stopcodon,$partial_start,$partial_stop) = &findCodons($atree,
									      $seqname,
									      $fmin,
									      $fmax,
									      $orient,$feat_name);


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
	    

	}
	else{
	    if($startcodon == -1){
		$feat_attrs->{$feat_name}->{'CS4'}++; #invalid start
	    }
	    else{
		$feat_attrs->{$feat_name}->{'CS3'}++;
	    }
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
	}
	else{
	    if($stopcodon == -1){
		$feat_attrs->{$feat_name}->{'CE4'}++; #invalid stop
	    }
	    else{
		$feat_attrs->{$feat_name}->{'CE3'}++;
	    }
	}
	if(exists $featstarts->{$feat_name} && $featstops->{$feat_name}){
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
	    $cluster_attrs->{'CS3'}++;
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
	    print "#Class CE3. Unaligned stops\n" if($debug);;
	    $cluster_attrs->{'CE1'}++;
	}
    }
    else{
	if($alignedstopcount == scalar(keys %$genes)){
	    #there is one annotated stop codon for each genome, but not all genomes use the same stop
	    print "#Class CE2. Inconsistent stops\n" if($debug);;
	    $cluster_attrs->{'CE2'}++;
	}
	else{
	    print "#Class CE3. Unaligned stops\n" if($debug);;
	    $cluster_attrs->{'CE3'}++;
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
    my $seqobj = $db->get_Seq_by_id($seqname);
    die "Bad coordinates $fmin-$fmax @_" if($fmin >= $fmax);
    my $codon_table = Bio::Tools::CodonTable->new(-id=>11);
    if($seqobj){
	if($orient eq '+'){
	    die "Bad coordinates $fmax extends past end of sequence" if($fmax >= $seqobj->length());
            my $seqlen = ($fmax-$fmin);
	    my $newobjs = $seqobj->trunc($fmin+1,$fmax);

	    my $encoding = 'C'x$newobjs->length();
	    my $newobj = new Bio::Seq::EncodedSeq(-seq=>$newobjs->seq(),
						  -encoding=>$encoding);
	    die if($newobj->length() != $seqlen);
	    if($codon_table->is_start_codon($newobj->subseq(1,3)) && ($codon_table->is_ter_codon($newobj->subseq($seqlen-3+1,$seqlen)))){
		my $protein_seq_obj = $newobj->translate(-codontable_id =>11);
		if($protein_seq_obj->length() == $seqlen/3){
		    return 1;
		}
		else{
		    print "#Unexpected sequence length ",$protein_seq_obj->length()," expecting ",$seqlen/3," from ORF $seqname $fmin-$fmax $orient\n" if($verbose);
		}
	    }
	    else{
		print "Possible alternative ORF on $seqname $fmin-$fmax,$orient has invalid start:",$newobj->subseq(1,3)," ",$codon_table->is_start_codon($newobj->subseq(1,3))," or stop:",$newobj->subseq($seqlen-3+1,$seqlen)," ",$codon_table->is_ter_codon($newobj->subseq($seqlen-3+1,$seqlen)),"\n" if($verbose);
	    }
	}
	else{
	    die if($orient ne '-');
            my $seqlen = ($fmax-$fmin);
	    my $newobj = $seqobj->trunc($fmin+1,$fmax);
	    die if($newobj->length() != $seqlen);
	    $newobj = $newobj->revcom();
	    #Check if valid start codon
	    if($codon_table->is_start_codon($newobj->subseq(1,3)) && ($codon_table->is_ter_codon($newobj->subseq($seqlen-3+1,$seqlen)))){
		my $protein_seq_obj = $newobj->translate(-codontable_id =>11);
		
		if($protein_seq_obj->length() == $seqlen/3){
		    return 1;
		}
		else{
		    print "#Unexpected sequence length ",$protein_seq_obj->length()," expecting ",$seqlen/3," from ORF $seqname $fmin-$fmax $orient\n" if($verbose);
		}
	    }
	    else{
		print "Possible alternative ORF on $seqname $fmin-$fmax,$orient has invalid start:",$newobj->subseq(1,3)," ",$codon_table->is_start_codon($newobj->subseq(1,3))," or stop:",$newobj->subseq($seqlen-3+1,$seqlen)," ",$codon_table->is_ter_codon($newobj->subseq($seqlen-3+1,$seqlen)),"\n" if($verbose);
	    }
	}
    }
    return 0;
}

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
		$newobj = $newobj->revcom();
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
	    $newobj = $newobj->revcom();
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
#
#Print members and attributes for a cluster
#$query is the longest member of a cluster
#Supported attributes
sub reportCluster{
    my($query,$mappedorgs,$mappedgenes,$unmappedgenes,$feat_attrs,$cluster_attrs,$seq_attrs,$new_orfs) = @_;
    if(scalar(keys %$mappedgenes)>0){
	if($COGoutputformat){
	    if(scalar(keys %$mappedgenes)>0){
		print "COG = $cluster_id, size ",scalar(keys %$mappedgenes), ", connections = 0, perfect = 0;\n";
		print "\t$features->{$query}->[5]\n";
		foreach my $organism (sort {$a cmp $b} keys %$mappedorgs){
		    foreach my $gene (sort {$features->{$a}->[1] <=> $features->{$b}->[1]} (keys %{$mappedorgs->{$organism}->{'features'}})){
			if($gene ne $query){
			    print "\t$features->{$gene}->[5]\n";
			}
		    }
		}
	    }
	}
	else{
	    my $classesstr = join(';',sort {$a cmp $b} keys %{$cluster_attrs});
	    print ">CLUSTER_$cluster_id num_seqs=",scalar(keys %$mappedorgs)," num_genes=",scalar(keys %$mappedgenes);
	    if(exists $mappedgenes->{$query}->{'alignments'}){
		print " num_alignments=",scalar(@{$mappedgenes->{$query}->{'alignments'}})," classes=$classesstr query=$query alignments=",join(',',@{$mappedgenes->{$query}->{'alignments'}});
	    }
	    print "\n";
	    print CFILE ">CLUSTER_$cluster_id num_seqs=",scalar(keys %$mappedorgs)," num_genes=",scalar(keys %$mappedgenes), " classes=$classesstr query=$query\n";

	    my $qfmin = $features->{$query}->[1];
	    my $qfmax = $features->{$query}->[2];
	    my $qseqname = $features->{$query}->[0];
	    my @mappedfeats;
	    foreach my $organism (sort {$a cmp $b} keys %$mappedorgs){
		my($start,$end) = &getspan($mappedgenes,keys %{$mappedorgs->{$organism}->{'features'}});
		my @ogenes = sort {$features->{$a}->[1] <=> $features->{$b}->[1]} (keys %{$mappedorgs->{$organism}->{'features'}});
		my @ocovs = map {sprintf("%.2f",$mappedgenes->{$_}->{'cov'}/$features->{$_}->[3])} (@ogenes); #%coverage over gene length
		my @oids  = map {sprintf("%.2f",$mappedgenes->{$_}->{'pid'}/$mappedgenes->{$_}->{'len'})} (@ogenes); #%id over aligned length
		my $classes;
		my $longestorf=0;
		my $longestpairc=0;
		foreach my $gene (@ogenes){
		    if(exists $feat_attrs->{$gene}){
		      foreach my $c (sort {$a cmp $b} keys %{$feat_attrs->{$gene}}){
			  #hack
			 if($c eq 'pairfreq'){
			     $classes->{'pairfreq='.$feat_attrs->{$gene}->{$c}}++;
			 }
			 else{
			     $classes->{$c}++;
			 }
		      }
                    }
		    $longestorf = ($features->{$gene}->[3] > $longestorf) ? $features->{$gene}->[3] : $longestorf;
		    $longestpairc = ($feat_attrs->{$gene}->{'pairfreq'} > $longestpairc) ? $feat_attrs->{$gene}->{'pairfreq'} : $longestpairc;
		}
		##
		#Report alternative start sites if they result in a longer ORF
		my @attrs = sort {$a cmp $b} keys %$classes;

		if(exists  $seq_attrs->{$organism}){
		    my $orfidx=0;
		    foreach my $alt (@{$seq_attrs->{$organism}}){
			#Report alternative starts or possible frameshifts
			if($alt =~ /alt_start/ || $alt =~ /alt_fs/){
			    #Only report if results in a longer ORF
			    my($astart,$aend) = ($alt =~ /alt_start=(\d+)-(\d+)/);
			    my($len) = ($alt =~ /len:(\d+)/);
			    my($orient) = ($alt =~ /orient:([\+\-])/);
			    my($apairfreq) = ($alt =~ /pairfreq:(\d+)/);
			    if(!$longer_altstarts || $len>=$longestorf){
				if(!$moreconsistent_altstarts || $apairfreq>=$longestpairc){
				    push @attrs,$alt;
				    push @mappedfeats,[[[$organism,$astart,$aend,$orient,($aend-$astart).'M']],"ALTgene$organism$orfidx",'gene'];
				    $orfidx++;
				    if($organism eq $qseqname){
					$qfmin = ($astart < $qfmin) ? $astart : $qfmin;	
					$qfmax = ($aend > $qfmax) ? $aend : $qfmax;
				    }
				}
				else{
				    print "#Skipping $alt codon freq $apairfreq<$longestpairc\n" if($debug);;
				}
			    }
			    else{
				print "#Skipping $alt len $len<$longestorf\n" if($debug);;
			    }
			}
			else{
			    print "#Unknown $alt" if($debug);;
			}
		    }
		}

		my @orients;
		my @names;
		foreach my $gene (@ogenes){
		    push @orients,"$features->{$gene}->[4]";
		    push @attrs,"aln_orient=$mappedgenes->{$gene}->{'relorient'}";
		    my $frame;
		    if($features->{$gene}->[4] eq '-'){
			$frame=($end%3)*-1;
		    }
		    else{
			$frame=$start%3;
		    }
		    push @attrs,"frame=$frame";
		    if(defined $features->{$gene}->[11]){
			push @names,"product=$features->{$gene}->[11]";
		    }
		}

		#Brief cluster output
		print CFILE join(',',@ogenes),
		"\tWGA$cluster_id",
		"\t$organism",
		"\tcov=",join(',',@ocovs),
		"\tpid=",join(',',@oids),
		"\tqcov=",sprintf("%.2f",$mappedorgs->{$organism}->{'qcov'}/($qfmax-$qfmin)),
		"\t$start-$end",
		"\t",join(',',@orients),
		"\t",$end-$start,
		"\t",join(',',@names),
		"\n";

		#Detailed output
		print join(',',@ogenes),
		"\tWGA$cluster_id",
		"\t$organism",
		"\tcov=",join(',',@ocovs),
		"\tpid=",join(',',@oids),
		"\tqcov=",sprintf("%.2f",$mappedorgs->{$organism}->{'qcov'}/($qfmax-$qfmin)),
		"\t$start-$end",
		"\t",join(',',@orients),
		"\t",$end-$start,
		"\t",join(';',@attrs,@names),
		"\n";
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
	    if($printalignments){
		print "#Printing query $query $qseqname,$qfmin,$qfmax\n" if($debug);
		my @isect = $atree->map($qseqname,$qfmin,$qfmax,"alignment");
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
		foreach my $align_name (@{$mappedgenes->{$query}->{'alignments'}}){
		    my $alni = $atree->getAlignedInterval($align_name,$feat2organism->{$query}) if($debug);
		    print "#QRYALN $align_name $alni->[1]\n" if($debug);
		    push @qryalns,[$align_name,$alni->[1]];
		}
		my @qryalns_sorted = sort {$a->[1] <=> $b->[1]} @qryalns;
		foreach my $a (sort {$a->[1] <=> $b->[1]} @qryalns){
		    my($align_name) = @$a;
		    #Check that new range is still within $alignment
		    print "#Checking the $qseqname,$qfmin,$qfmax,$align_name is within range\n" if($debug);;
		    my @isect = $atree->intersect($qseqname,$qfmin,$qfmax,$align_name);
		    my $printfmin;
		    my $printfmax;
		    foreach my $aln (@isect){
			if($aln->[1] eq $qseqname && $aln->[0] eq $align_name){
			    #print join(',',@$aln),"\n";
			    $printfmin = $aln->[2];
			    $printfmax = $aln->[3];
			    print "#Resetting print range to $printfmin-$printfmax from $qfmin-$qfmax\n" if($debug);
			}
		    }
		    if(defined $printfmin && defined $printfmax){
			print "CLUSTER_$cluster_id ALIGNMENT:$align_name\n";
			my($colstart,$colend) = AlignmentTree::coordstocolumn($atree->{_alignments}->{$align_name}->[0],$qseqname,$printfmin,$printfmax);
			$atree->printAlignment($align_name,$colstart,$colend,$db,\@mappedfeats);
		    }
		    else{
			die;
		    }
		}
	    }
	    print "\n";
	}
    }
    else{
	#No genes in cluster
	die;
    }
}

sub findSingletons{
    my($atree,$mapped,$unmapped,$subsumed) = @_;
    my $singletons = {};
    foreach my $feat_name (keys %$features){
	my $fmin = $features->{$feat_name}->[1];
	my $fmax = $features->{$feat_name}->[2];
	if(! exists $mapped->{$feat_name}){
	    die if(exists $mapped->{$feat_name});
	    my $classes = &annotateSingletons($atree,$features->{$feat_name}->[0],$feat_name,$fmin,$fmax);
	    if(exists $unmapped->{$feat_name}){
		my $query=$feat_name;
		my($mappedorgs,$mappedgenes,$unmappedorgs,$unmappedgenes) = &buildCluster($atree,$query);
		my($feat_attrs,$cluster_attrs,$codons) = &annotateCluster($atree,$mappedgenes,$mappedorgs);
		my $new_orfs = &findnewORFs($db,$atree,$mappedorgs,$mappedgenes,$codons);
		if(scalar(keys %$new_orfs)){
		    my $seq_attrs = {};	 
		    &reportCluster($query,$mappedorgs,$mappedgenes,$unmappedgenes,$feat_attrs,$cluster_attrs,$seq_attrs,$new_orfs);
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
		if($COGoutputformat){}
		else{
		    print "#SINGLETON $feat_name len:$features->{$feat_name}->[3]\tbest_cluster:$unmapped->{$feat_name}->{'WGA_cluster'}\tcov:";
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
		    print "\n";
		}
	    }
	    else{
		if(exists $subsumed->{$feat_name}){
		    print "#DELETED $feat_name\n";
		}
		else{
		    if($COGoutputformat){}
		    else{
			print "#SINGLETON $feat_name len:$features->{$feat_name}->[3] ",join(' ',@$classes);
			if(defined $features->{$feat_name}->[11]){
			    print " product=$features->{$feat_name}->[11]";
			}
			print "\n";
		    }
		    $nohit++;
		    $singletons->{$feat_name}++;
		}
	    }
	}
	else{
	    #Mapped ORF, not a singleton
	}
    }

    return $singletons;
}


###############################
#General utility funcs
sub getspan{
    my($features) = shift;
    my @coords;
    foreach my $gene (@_){
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
	if($ei){
	    ($stop_s,$stop_e) = AlignmentTree::columntocoords($ei,$stopcol,$stopcol+2);
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
    my $seqobj = $db->get_Seq_by_id($seqname);
    if(!$seqobj){
	print "Can't find $seqname\n";
	return;
    }
    my $startcodon=undef;
    my $stopcodon=undef;
    my $is_partial_start=0;
    my $is_partial_stop=0;
    my $aln_orient=undef;
    if($orient eq '+'){
	if(!$codon_table->is_start_codon($seqobj->subseq($fmin+1,$fmin+2+1))){ #bioperl is 1-base coordinates
	    print "#Bad start codon $fname,$seqname,$fmin,$fmax,$orient codon $fmin+1,$fmin+2+1 ",$seqobj->subseq($fmin+1,$fmin+2+1)," aln_orient:$aln_orient\n" if($verbose || $debug);
	    return -1;
	} 
	else{
	    #Find start codon + strand
	    $startcodon = &getAlignedCols($atree,$seqname,$fmin,$fmin+3);
	}
	
	if(!$codon_table->is_ter_codon($seqobj->subseq($fmax-3+1,$fmax))){
	    print "#Bad stop $fname,$seqname,$fmin,$fmax,$orient codon $fmax-3+1,$fmax ",$seqobj->subseq($fmax-3+1,$fmax)," aln_orient:$aln_orient\n" if($verbose || $debug);
	    return -1;
	}
	else{
	    #Find stop codon - strand
	    $stopcodon = &getAlignedCols($atree,$seqname,$fmax-3,$fmax);
	}
	#Check if in pmark spacer adjacent to contig boundary
	my $startregion = $seqobj->subseq($fmin-length($PMARK_SPACER),$fmin+length($PMARK_SPACER));
	my $stopregion = $seqobj->subseq($fmax-length($PMARK_SPACER),$fmax+length($PMARK_SPACER));
	if($startregion =~ /$PMARK_SPACER/){
	    $is_partial_start=1;
	}
	if($stopregion =~ /$PMARK_SPACER/){
	    $is_partial_stop=1;
	}
	
    }
    else{
	die "Bad orient $orient" if($orient ne '-');
	if(!$codon_table->is_start_codon(revcom($seqobj->subseq($fmax-3+1,$fmax))->seq())){
	    print "#Bad start codon $fname,$seqname,$fmin,$fmax,$orient codon $fmax-3+1,$fmax ",revcom($seqobj->subseq($fmax-3+1,$fmax))->seq()," aln_orient:$aln_orient\n" if($verbose || $debug);
	    return -1;
	} 
	else{
	    #Find start codon on - strand
	    $startcodon = &getAlignedCols($atree,$seqname,$fmax-3,$fmax);
	}
	if(!$codon_table->is_ter_codon(revcom($seqobj->subseq($fmin+1,$fmin+3))->seq())){
	    print "#Bad stop codon $fname,$seqname,$fmin,$fmax,$orient codon $fmin+1,$fmin+3 ",revcom($seqobj->subseq($fmin+1,$fmin+3))->seq()," aln_orient:$aln_orient\n" if($verbose || $debug);
	    return -1;
	} 
	else{
	    #Find stop codon on - strand
	    $stopcodon = &getAlignedCols($atree,$seqname,$fmin,$fmin+3);
	}
        #Check if in pmark spacer adjacent to contig boundary
	my $stopregion = $seqobj->subseq($fmin-length($PMARK_SPACER),$fmin+length($PMARK_SPACER));
	my $startregion = $seqobj->subseq($fmax-length($PMARK_SPACER),$fmax+length($PMARK_SPACER));
	if($startregion =~ /$PMARK_SPACER/){
	    $is_partial_start=1;
	}
	if($stopregion =~ /$PMARK_SPACER/){
	    $is_partial_stop=1;
	}
    }
    return ($startcodon,$stopcodon,$is_partial_start,$is_partial_stop);
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
    my ($mmatrix,$seqmatrix,$names) = $atree->getAlignmentMatrix($aln,$startcol,$endcol,$db);
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
    foreach my $r (keys %$results){
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
	my $alni = &getAlignment($atree,$aln,$seq);
	my($fsstart,$fsend) = AlignmentTree::columntocoords($alni,$reloffset,$reloffset);
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
		    print "#ALT col:$reloffset coord:$fsstart-$fsend base:$refchar freq:$freqchar->{$refchar} $seq:$qrychar $freqchar->{$qrychar} fstype:$fstype\n" if($debug);
		    push @edits,[$fsstart,$refchar,$qrychar,$reloffset,$fstype];
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
		    print "#WARNING Overlapping aligned region found for $seqname,$fmin,$fmax. $align_name and $ret->[1]\n" if($verbose);
		}
		my @res= AlignmentTree::coordstocolumn($alni,$seqname,$fmin,$fmax);
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
			    my $seqobj = $db->get_Seq_by_id($seqname);
			    if($seqobj){
				die "Can't find sequence $seqname obj:$seqobj" if(!defined $seqobj);
				my ($neworf,$callorient) = &callORF($seqobj,$start,$end,$orient);
				if(length($neworf)>$MINORF){
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
				    print "Skipping short ORF ",length($neworf)," <$MINORF $start,$end,$orient\n" if($debug);
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
    my($atree,$seqname,$neighborseqs,$startcodon,$stopcodon) = @_;
    #Short circuit for testing, return any neighbor
    return @$neighborseqs[0];
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
	
	my $fsvars;
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
		my ($startcol,$stopcol) = AlignmentTree::coordstocolumn($atree->{_alignments}->{$aln->[0]}->[0],$seqname,$aln->[2],$aln->[3]);
		my $sv = &reportVariants($atree,$db,$aln->[0],$seqname,$startcol,$stopcol,$nearestseq);
		foreach my $v (@$sv){
		    if(abs($netfs) > $FS_THRESHOLD){
			#Short circuit
		    return undef;
		}
		    else{
			if($v->[4] != 0){
			    print "#FSVAR $seqname ",join(',',@$v),"\n" if($debug);
			    push @$fsvars,$v;
			    $netfs += $v->[4];
			}
		    }
		}
	    }
	}
	return ($fsvars,$netfs);
    }
}
##############################
#Report annotations that can be reconciled with a frameshift
#These may indicate a sequencing error or an authentic frameshift
#
#Method:

#Consider all clusters with inconsistent start and/or stops and report
#ORFs that are consistent after a frameshift. Any indels occuring
#within a distance d upstream from an annotated stop are considered in
#series as a possible frameshift.  Uses currently annotated start
#codons and also any alternative start sites found with checkStarts()

#Consider all start codons. Find any indels that would shift frame 
#or cause a premature stop (PM that result in a stop, indels of len != mod 3)
# 


#Reports 
#(1) frameshift location and indel
#(2) resulting ORF,length, and translation
#(3) any overlapping genes that would be subsumed by the new ORF

#fs is signed locations of the frameshift relative to the sequence start
#eg. +10 is a forward frameshift 10 bp downstream from translation start
#    -9 is a backward frameshift 9 bp downstream from translation start

sub checkFrameshifts_new{
    my ($db,$codons,$seqs,$seq_attrs) = @_;
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
    
    foreach my $seq (keys %$seqs){ 
	#Consider all start codons and find corresponding stop
	if(scalar(keys %{$seqs->{$seq}->{'features'}})>1){
	    foreach my $pcodon (keys %{$codons->{'pairs'}}){
		my($startcodon,$stopcodon) = split(/:/,$pcodon);
		#If using this stop codon, look for a frameshift
		if(exists $codons->{'stops'}->{$seq}->{$stopcodon}){
		print "Looking for FS on $seq using $startcodon , $stopcodon\n";
		my $fswin_min;
		my $fswin_max;
		my $start_s;
		my $start_e;
		my $stop_s;
		my $stop_e;
		#$codon is a tuple of alignment,aligned_column
		my($startcol,$aln_s) = split(/$CODON_DELIM_REGEX/,$startcodon);
		#find corresponding stop
		my($stopcol,$aln_e) = split(/$CODON_DELIM_REGEX/,$stopcodon);
		
		my $si = &getAlignment($atree,$aln_s,$seq);
		my $ei = &getAlignment($atree,$aln_e,$seq);
		if($si){
		    ($start_s,$start_e) = AlignmentTree::columntocoords($si,$startcol,$startcol+2);
		    if($ei){
			($stop_s,$stop_e) = AlignmentTree::columntocoords($ei,$stopcol,$stopcol+2);
			$fswin_max=$stop_e;
			$fswin_min=$stop_e-$FS_WINDOW;
			my $origorflen = abs($stop_e-$start_s);
			if($atree->contains($aln_e,$seq,$fswin_min,$fswin_max)){
			    print "Contained within alignment\n";
			    my $alni = &getAlignment($atree,$aln_e,$seq);
			    my @res= AlignmentTree::coordstocolumn($atree->{_alignments}->{$aln_e}->[0],$seq,$fswin_min,$fswin_max);
			    my $startcol = $res[0];
			    my $stopcol = $res[1];
			    my $indels = &reportVariants($atree,$db,$aln_e,$seq,$startcol,$stopcol);
			    #print "#$codon freq:$codons->{'starts'}->{$seq}->{$codon} orient:$orient ORF:$orfstart-$orfend len:",$orfend-$orfstart," aalen:",length($origorf),"\n" if($debug);;
			    foreach my $indel (@$indels){
				my($fsstart,$indelbase,$origbase) = @$indel;
				#die "$fsstart < $start-$end\n" if($fsstart < $start);
				my $fsloc;
				my $fstype;
				if($indelbase eq '-'){
				    #Reverse frameshift
				    $fsloc = ($fsstart-$start_s) * -1; #
				    $fstype = "F";
				}
				elsif($origbase eq '-'){
				    #Forward frameshift
				    $fsloc = ($fsstart-$start_s);
				    $fstype = "R";
				}
				elsif($origbase ne '-' && $indelbase ne '-'){
				    #point mutation
				    #Not handled
				}
				else{
				    die "Invalid combination $fsstart,$indelbase,$origbase";
				}
				$fsstart = ($fsstart < $start_s) ? $start_s : $fsstart;
				#if($origbase ne 'X'){
				#   my $base = substr($seqobj->seq(),$fsstart,1);
				#    die if($base ne $origbase);
				#}
				print "#Edit @ $fsstart $origbase->$indelbase\n";# if($debug);;
				print "#Attempting frameshift for ORF $start_s-$start_e @ loc $fsloc\n";
				my $orient = $alni->[3]; #TODO check this
				my $seqobj = $db->get_Seq_by_id($seq);
				
				
				my ($neworf,$neworient) = &callORF($seqobj,$start_s,$start_e,$orient,[$fsloc]);
				print "Length neworf:",length($neworf)*3," vs. orig annotation ",$origorflen,"\n";
				#Only consider frameshifts which extend the ORF
				if((length($neworf)*3) > $origorflen){
				    print "#NEWORF ",length($neworf)*3," vs. $origorflen\n";
				    my $fmin;
				    my $fmax;
				    
				    #Get start coordinate and codon location
				    
				    #Get stop coordinate and codon location
				    if($neworient eq '+'){
					$fmin=$start_s;
					$fmax=$start_s+(length($neworf)*3);
					if($fstype eq 'F'){ #is forward frameshift
					    #$fmax -= 1;
					}
					else{
					    #$fmax += 1;
					}
				    }
				    elsif($neworient eq '-'){
					$fmin=$start_e-(length($neworf)*3);
					$fmax=$start_e;
					if($fstype eq 'R'){ #is forward frameshift
					    #$fmin -= 1;
					}
					else{
					    #$fmin += 1;
					}
				    }
				    print "Finding codons for ORF $fmin-$fmax,$orient $fstype vs. orig $start_s,$start_e - $stop_s,$stop_e\n";
				    my($new_s,$new_e) = &findCodons($atree,
								    $seq,
								    $fmin,
								    $fmax,
								    $neworient);
				    my($mcol,$align_name) = (@$new_s);
				    my $new_startcodon = $mcol.$CODON_DELIM.$align_name;
				    ($mcol,$align_name) = (@$new_e);
				    my $new_stopcodon = $mcol.$CODON_DELIM.$align_name;
				    print "Codons: $new_startcodon $new_stopcodon on $seq\n";
				    
				    die "Orig: $pcodon New:$new_startcodon:$new_stopcodon" if($new_startcodon != $startcodon);
				    
				    #Add this ORF to list of alternatives for this genome
				    if(! exists $codons->{'starts'}->{$seq}->{$startcodon}){
					$codons->{'starts'}->{$seq}->{$startcodon}++;
				    }
				    if(! exists $codons->{'stop'}->{$seq}->{$new_stopcodon}){
					$codons->{'stop'}->{$seq}->{$new_stopcodon}++;
				    }
				    if(! exists $codons->{'pairs'}->{$startcodon.':'.$new_stopcodon}->{'orgs'}->{$seq}){
					$codons->{'pairs'}->{$startcodon.':'.$new_stopcodon}->{'orgs'}->{$seq} = [$fmin,$fmax,$neworient,0,$fsstart];
				    }
				    else{
					print "#Codon pair $startcodon.':'.$new_stopcodon already exists on $seq ",join(',',@{$codons->{'pairs'}->{$startcodon.':'.$new_stopcodon}->{'orgs'}->{$seq}}),"\n";
				    }
				}
			    }
			}
			else{
			    print "Skipping frameshift check. $fswin_min - $fswin_max is not contained within alignment $aln_e\n" if($debug);
			}
		    }
		}
		}
	    }
	}
    }
}
    
sub checkFrameshifts{
    my ($db,$seqs,$genes,$codons,$seq_attrs) = @_;

    #Save list of all start codons $codon->$freq
    my $starts = {};
    my $stops = {};
    foreach my $seq (keys %{$codons->{'starts'}}){
	foreach my $codon (keys %{$codons->{'starts'}->{$seq}}){
	    $starts->{$codon} += $codons->{'starts'}->{$seq}->{$codon};
	}
    }
    
    foreach my $seq (keys %$seqs){ #sort on span shortest to longest
	my $neworflist = {}; #candidate new orfs
	#Consider all start codons
	#Find location of indels on or upstream from the stop (< bp cutoff)
	foreach my $codon (keys %$starts){
	    #check if $codon is aligned and upstream
	    #$codon is a tuple of alignment,aligned_column
	    my($startcol,$aln) = split(/$CODON_DELIM_REGEX/,$codon);
	    my $gapped=1;
	    #check is $col,$col+3 is gapped, return start coordinate on the genome
	    my $i = &getAlignment($atree,$aln,$seq);
	    if($i){
		die "Cannot find alignment $aln that contains $seq" if(!$i);
		my($start,$end) = AlignmentTree::columntocoords($i,$startcol,$startcol+2);
		$gapped = ($end-$start == 3) ? 0 : 1;		
		if(!$gapped){ #codon is aligned, attempt to call and ORF with a frameshift or PM
		    my @neworfs;
		    my $orient = $i->[3];
		    my $seqobj = $db->get_Seq_by_id($seq);
		    my ($origorf,$origorforient) = &callORF($seqobj,$start,$end,$orient);
		    #Check if there is an ORF
		    if(defined $origorf){
			my ($orfstart,$orfend,$orflen);
			if($origorforient eq '+'){
			    $orfstart=$start;
			    $orfend=$start+(length($origorf)*3);
			}
			else{
			    $orfstart=$end-(length($origorf)*3);
			    $orfend=$end;
			}
			$orflen = $orfend-$orfstart;
			print "#Looking for stop codon $orfend-3,$orfend in $aln for ORF $orfstart-$orfend\n#$origorf\n" if($debug);;
			#Check if alignment contains both current start and stop codon
			if($atree->contains($aln,$seq,$orfstart,$orfstart+3) &&
			   $atree->contains($aln,$seq,$orfend-3,$orfend)){
			    my $stopcol;
			    if($origorforient eq '+'){
				my @res= AlignmentTree::coordstocolumn($atree->{_alignments}->{$aln}->[0],$seq,$orfend-3,$orfend);
				$stopcol = $res[0];
			    }
			    else{
				my @res= AlignmentTree::coordstocolumn($atree->{_alignments}->{$aln}->[0],$seq,$orfstart,$orfstart+3);
				$stopcol = $res[0];
			    }
			    my $indels = &reportVariants($atree,$db,$aln,$seq,$startcol,$stopcol);
			    print "#$codon freq:$codons->{'starts'}->{$seq}->{$codon} orient:$orient ORF:$orfstart-$orfend len:",$orfend-$orfstart," aalen:",length($origorf),"\n" if($debug);;
			    foreach my $indel (@$indels){
				my($fsstart,$indelbase,$origbase) = @$indel;
				#die "$fsstart < $start-$end\n" if($fsstart < $start);
				$fsstart = ($fsstart < $start) ? $start : $fsstart;
				#if($origbase ne 'X'){
				#   my $base = substr($seqobj->seq(),$fsstart,1);
				#    die if($base ne $origbase);
				#}
				print "#Edit @ $fsstart $origbase->$indelbase\n" if($debug);;
				my($neworf,$newseq,$neworfstart,$neworfend,$neworflen,$neworforient);
				if($orient eq '+'){
				    #build new string on leading strand
				    my $seqlen = ($seqobj->length()>$MAXORFLEN) ? $start+$MAXORFLEN : $seqobj->length();
				    my $newobj = $seqobj->trunc($start+1,$seqlen);
				    die "bad trunc" if($newobj->subseq(1,3)!=$seqobj->subseq($start+1,$start+3));
				    my $newcoord = $fsstart-$start;
				    die "$newcoord $fsstart $start" if($newcoord<0);
				    #my $mutation;

				    if($indelbase eq '-'){
					#remove a base frameshift
					$newseq = substr($newobj->seq,0,$newcoord-1);
					$newseq .= substr($newobj->seq(),$newcoord+1,$newobj->length());
				    }
				    elsif($origbase eq '-'){
					#add a base, frameshift
					$newseq = substr($newobj->seq,0,$newcoord);
					$newseq .= $indelbase.substr($newobj->seq(),$newcoord+1,$newobj->length());
				    }
				    else{
					#point mutation
					$newseq = substr($newobj->seq,0,$newcoord);
					$newseq .= $indelbase.substr($newobj->seq(),$newcoord+1,$newobj->length());
					die "$start $fsstart $newcoord target:$indelbase found:".substr($newseq,$newcoord,1)."orig:$origbase" if($indelbase ne substr($newseq,$newcoord,1));
				    }
				}
				else{
				    #build new string on leading strand
				    my $seqlen = ($end>$MAXORFLEN) ? $end-$MAXORFLEN : 1;
				    my $newobj = $seqobj->trunc(1,$end+1);
				    $newobj = $newobj->revcom();
				    my $newcoord = $end-$fsstart;
				    die "$fsstart $end" if($newcoord <= 0);
				    if($indelbase eq '-'){
					#remove a base frameshift
					$newseq = substr($newobj->seq,0,$newcoord-1);
					$newseq .= substr($newobj->seq(),$newcoord+1,$newobj->length());
				    }
				    elsif($origbase eq '-'){
					#add a base, frameshift
					$newseq = substr($newobj->seq,0,$newcoord);
					$newseq .= $indelbase.substr($newobj->seq(),$newcoord+1,$newobj->length());
				    }
				    else{
					#point mutation
					$newseq = substr($newobj->seq,0,$newcoord);
					$newseq .= $indelbase.substr($newobj->seq(),$newcoord+1,$newobj->length());
					die "$start $fsstart $newcoord target:$indelbase found:".substr($newseq,$newcoord,1)."orig:$origbase" if($indelbase ne substr($newseq,$newcoord,1));
				    }
				}
				my $newobj = Bio::Seq->new(-seq=>$newseq);
				($neworf,$neworforient) = &callORF($newobj,0,$newobj->length()-1,$orient);
				if(! defined $neworf){
				    #TODO bad start, skip till next start
				    last;
				}
				else{
				    if($neworforient eq '+'){
					$neworfstart = $start;
					$neworfend = $start + (length($neworf)*3);
				    }
				    else{
					$neworfstart=$end-(length($neworf)*3);
					$neworfend=$end; 
				    }
				    $neworflen = $neworfend-$neworfstart;
				}
				print "#Considering FS ORF $neworfstart-$neworfend:$neworflen $neworf ",substr($newseq,0,$neworflen),"\n#original:$orfstart-$orfend $orflen\n" if($debug);;
				#If ORF passes minimum length and has a new stop codon
				if(length($neworf)>$MINORF 
				   && $neworfend != $orfend
				   && $neworflen > $orflen){
				    print "#New orf on $seq len:",length($neworf)," ORF:$neworfstart-$neworfend\n" if($debug);;
				    if($atree->contains($aln,$seq,$neworfend,$neworfend+1)){
					my @res= AlignmentTree::coordstocolumn($atree->{_alignments}->{$aln}->[0],$seq,$neworfend,$neworfend+1);
					$stopcol = $res[0];
				    }
				    else{
					#new stop runs outside alignment
					$stopcol = '?';
				    }
				    $neworflist->{$seq}->{$neworfstart.$CODON_DELIM.$neworfend}->{'orf'} = $neworf;
				    $neworflist->{$seq}->{$neworfstart.$CODON_DELIM.$neworfend}->{'start'} = $neworfstart;
				    $neworflist->{$seq}->{$neworfstart.$CODON_DELIM.$neworfend}->{'end'} = $neworfend;
				    $neworflist->{$seq}->{$neworfstart.$CODON_DELIM.$neworfend}->{'orient'} = $neworforient;
				    $neworflist->{$seq}->{$neworfstart.$CODON_DELIM.$neworfend}->{'fsstart'} = $fsstart;
				    $neworflist->{$seq}->{$neworfstart.$CODON_DELIM.$neworfend}->{'orig'} = $origbase;
				    $neworflist->{$seq}->{$neworfstart.$CODON_DELIM.$neworfend}->{'edit'} = $indelbase;
				    $neworflist->{$seq}->{$neworfstart.$CODON_DELIM.$neworfend}->{'startcol'} = $startcol;
				    $neworflist->{$seq}->{$neworfstart.$CODON_DELIM.$neworfend}->{'stopcol'} = $stopcol;
				}
			    }
			}
		    }
		}
	    }
	}
	#Get list of new orfs
	foreach my $nf (keys %{$neworflist->{$seq}}){
	    my $neworf = $neworflist->{$seq}->{$nf}->{'orf'};
	    my $start = $neworflist->{$seq}->{$nf}->{'start'};
	    my $end = $neworflist->{$seq}->{$nf}->{'end'};
	    my $orient = $neworflist->{$seq}->{$nf}->{'orient'};
	    my $fsstart = $neworflist->{$seq}->{$nf}->{'fsstart'};
	    my $origbase = $neworflist->{$seq}->{$nf}->{'orig'};
	    my $indelbase = $neworflist->{$seq}->{$nf}->{'edit'};
	    #TODO is longer than annotated ORFs on this seq
	    #TODO is consistent with another annotated stop codon
	    if(!exists $seq_attrs->{$seq}){
		$seq_attrs->{$seq} = [];
	    }
	    push @{$seq_attrs->{$seq}},"alt_fs=$start-$end,orient:$orient,len:".($end-$start).",edit: $fsstart $origbase -> $indelbase startcol:$neworflist->{$seq}->{$nf}->{'startcol'},stopcol:$neworflist->{$seq}->{$nf}->{'stopcol'}";
	}
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


