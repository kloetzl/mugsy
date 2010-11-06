#!/usr/bin/perl
######################
=head1 NAME

mapfeatures - derives a set of mapped features according to a
multiple sequence alignment. Reports on the consistency of
annotated features in the mapping.

=head1 USAGE

mapfeatures.pl alignments.index seqs.fasta < features.txt 

(1) alignment.index - A index file that is generated with the
indexing scripts featureindex.pl,mafindex.pl,xmfaindex.pl. The index
should contain at least one set of aligned regions such as produced
by a whole genome aligner like MUGSY, Mauve, or TBA

(2) features.txt - A space delimited file consisting of 
feature_id sequence_id fmin fmax strand

=head1 SYNOPSIS
#############
#APPLICATIONS
#############

1)Reporting orthologs according to a whole genome alignment

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

=cut


use strict;
use lib '/usr/local/projects/angiuoli/developer/sangiuoli/mugsy/trunk/mapping';

use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

#Bioperl is used only for translation machinery
use Bio::Perl;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::Tools::CodonTable;
#use Bio::LiveSeq::Mutation; tried this but couldn't get to work properly
#Default cutoffs

use AlignmentTree;

my %options;
my $results = GetOptions (\%options, 
			  'coverage|c=s',
			  'query_coverage|q=s',
			  'identity|i=s',
			  'minorflen=s',
			  'maxorflen=s',
			  'cogformat=s',
			  'printalignments=s',
			  'skipaltstarts',
			  'skipneworfs',
			  #'skipframeshifts',
			  'debug|d=s') || pod2usage(-verbose => 1);

pod2usage(-verbose=>1) if($options{'help'});

#TODO make configurable
my $coverage_cutoff = (exists $options{'coverage'}) ?  $options{'coverage'} : 0.7;
my $query_coverage_cutoff = (exists $options{'query_coverage'}) ?  $options{'query_coverage'} : 0;
my $pid_cutoff= (exists $options{'identity'}) ?  $options{'identity'} : 0.6;
print STDERR "Using coverage cutoff:$coverage_cutoff identity:$pid_cutoff query_coverage:$query_coverage_cutoff\n";

my $MAXORFLEN = (exists $options{'maxorflen'}) ? $options{'maxorflen'} : 30000;
#Flag for checking consistent start,stop
#Assumes input features are genes
my $doconsistencychecks=1;
#Report new ORFs using aligned start codons
my $dofindneworfs = 1;
my $autocorrect=0;
my $MINORF= $options{'minorflen'} || 50;#aa
#my $autofixunmapped=0; #TODO, use consistency checks to fix annotations of unmapped genes

#Only report alternative start codons that
#results in a longer ORF
my $longer_altstarts=1;
#Only report alternative start codons that
#appear more frequently in the aligned genoems
my $freq_altstarts=1;
my $freq_altstops=0;

my $aligntoken="WGA";

#Output flags
my $COGoutputformat=(exists $options{'cogformat'}) ? $options{'cogformat'} : 0;
my $printskipped=1;
my $printalignments=(exists $options{'printalignments'}) ? $options{'printalignments'} : 0;

#Debugging flags
my $checkbadlen=0;
my $debug=$options{'debug'};
my $verbose=0;

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
while(my $line=<STDIN>){
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
		foreach my $feat_name (keys %$mappedgenes){
		    if($mappedgenes->{$feat_name}->{'morient'}!=0){
			print "#Class O1. Mis-matched orientation on $feat_name\n" if($debug);;
			$feat_attrs->{$feat_name}->{'classO1'} = $mappedgenes->{$feat_name}->{'morient'};
		    }
		}
		
		#A consistent cluster has both consistent start and stop codons (labeled CS1,CE1)
		my $consistent=(exists $cluster_attrs->{'CS1'} && exists $cluster_attrs->{'CE1'});
		
		if(!$consistent){
		    #(2) Attempt to adjust start sites
		    #Using aligned codons as possible alternative start sites
		    #Print out alternative (fixed) gene coordinates		    
		    print "#Looking for new starts\n" if($debug);;
		    if(!defined $options{'skipaltstarts'}){
			&checkStarts($db,$codons,[keys %$mappedorgs],$seq_attrs);
		    }
		    #(3) Attempt to resolve inconsistencies by frameshifting
		    #Using aligned codons as start sites and indels adjacent to stops
		    #as possible locations of the frameshift
		    #Long ORFs may indicate a sequencing error or authentic frameshift
		    #Reports 
		    #-frameshift location and indel
		    #-new ORF,translation
		    #-any deleted genes
		    #if(!defined $options{'skipframeshifts'}){
		    #&checkFrameshifts($db,$mappedorgs,$mappedgenes,$codons,$seq_attrs);
		    #}
		}
	    }
	
	    #We have a good cluster, save it
	    #Save the cov,pid in master list of mapped genes
	    foreach my $feat_name (keys %$mappedgenes){
		die "Feature $feat_name already mapped" if(exists $mapped->{$feat_name});
		$mapped->{$feat_name}->{'cov'}=$mappedgenes->{$feat_name}->{'cov'}/$features->{$feat_name}->[3];
		$mapped->{$feat_name}->{'pid'}=$mappedgenes->{$feat_name}->{'pid'}/$mappedgenes->{$feat_name}->{'len'};
		delete $unmapped->{$feat_name};
	    }

	    #Print cluster
	    &reportCluster($query,$mappedorgs,$mappedgenes,$unmappedgenes,$feat_attrs,$cluster_attrs,$seq_attrs,$new_orfs);
	    my $classesstr = join(';',sort {$a cmp $b} keys %{$cluster_attrs});
	    $classes_sum->{$classesstr}->{'ngenes'} +=scalar(keys %$mappedgenes);
	    $classes_sum->{$classesstr}->{'nclusters'}++;

	    $validcluster++;
	    if($COGoutputformat){}
	    else{
		print "#VALID\tWGA$cluster_id\tNum_organisms=",scalar(keys %$mappedorgs)+1,
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
	    foreach my $gene (@ogenes){
		if(exists $feat_attrs->{$gene}){
		    foreach my $c (sort {$a cmp $b} keys %{$feat_attrs->{$gene}}){
			$classes->{$c}++;
		    }
		}
		$longestorf = ($features->{$gene}->[3] > $longestorf) ? $features->{$gene}->[3] : $longestorf;
	    }
	    ##

	    my @attrs = sort {$a cmp $b} keys %$classes;
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
			print "#alt $astart,$aend,$aorient,$alen\n" if($debug);;
			die "$alt" if(!$astart || !$aend || !$aorient || !$alen);
			if(!$longer_altstarts || $alen>$longestorf){
			    push @alts,["ALTSTARTgene$organism$orfidx",$astart,$aend,$aend-$astart,$aorient];
			    $orfidx++;
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
print "CM1 - multiple fragments spanned\n";
print "CO1 - mismatch orientation\n";
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
    
    #List of organism_ids in the current cluster 
    my $mappedorgs = {}; #passes cutoffs
    my $unmappedorgs = {}; #do not pass cutoffs
    
    #List of genes in the current cluster
    my $mappedgenes = {}; #passes cutoffs
    my $unmappedgenes = {}; #do not pass cutoffs
    
    #Contains list of annotations that are overlapping in an alignment
    my $alnfeats = {};
    my $alnorgs = {};

    #Map of feat_name->organism_name
    my $feat2organism = {};
    
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
	}
    }
    
    if($COGoutputformat){}
    else{
	print "#QUERY=$query coords=$qfmin-$qfmax len=$features->{$query}->[3] strand=$features->{$query}->[4] Num_alignments=",scalar(keys %$goodalignments),"\n";
    }
    
    foreach my $r (sort {   #Sort alphanumeric on alignment_name, secondary on align_start
	if($a->[0] eq $b->[0]){
	    $a->[2] <=> $b->[2];
	}
	else{
	    $a->[0] cmp $b->[0];
	}
    } @isect){
	my $feat_name = $r->[0];
	my $seqname = $r->[1];
	my $align_name = $r->[5];
	#Check if we want to consider this alignment
	if(exists $goodalignments->{$align_name}){
	    my($alnobj,$bv,$width) = $atree->getAlignment($align_name);
	    $feat_name =~ s/gene\://;
	    if(!exists $features->{$feat_name}){
		print STDERR "#Bad feature found $feat_name. Not in input file. Skipping\n";
		next;
	    }
	    #Capture some stats on the matching genes
	    #TODO the cov,pid stats assume non-overlapping alignments
	    if($query ne $feat_name){
		#Only report genes that have not been mapped
		if(!exists $mapped->{$feat_name} && !exists $deleted->{$feat_name}){
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
	    $unmappedorgs->{$feat2organism->{$feat_name}}->{'features'}->{$feat_name}++;
	    $unmappedorgs->{$feat2organism->{$feat_name}}->{'features'}->{$feat_name}++;
	    $unmappedorgs->{$feat2organism->{$feat_name}}->{'qcov'} = $alnorgs->{$feat2organism->{$feat_name}}->{'qcov'};
	    
	    $unmappedgenes->{$feat_name}->{'cov'} = $alnfeats->{$feat_name}->{'cov'};
	    $unmappedgenes->{$feat_name}->{'fmin'} = $alnfeats->{$feat_name}->{'fmin'};
	    $unmappedgenes->{$feat_name}->{'fmax'} = $alnfeats->{$feat_name}->{'fmax'};
	    $unmappedgenes->{$feat_name}->{'pid'} = $alnfeats->{$feat_name}->{'pid'};
	    $unmappedgenes->{$feat_name}->{'len'} = $alnfeats->{$feat_name}->{'len'};
	    $unmappedgenes->{$feat_name}->{'relorient'} = $alnfeats->{$feat_name}->{'relorient'};

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
	my($startcodon,$stopcodon) = &findCodons($atree,
						 $seqname,
						 $fmin,
						 $fmax,
						 $orient,$feat_name);

	if(ref $startcodon){
	    my($mcol,$align_name) = (@$startcodon);
	    my $token = $mcol.'?'.$align_name;
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
	    my $token = $mcol.'?'.$align_name;
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
    }

    if(scalar(keys %$starts)==1){
      my @start = keys %$starts; 
      if($starts->{$start[0]}==scalar(keys %$genes)){
	print "#Class CS1. Consistent starts\n" if($debug);;
	$cluster_attrs->{'CS1'}++;
      }
      else{
        print "#Class CS3. Unaligned starts ",$starts->{$start[0]}, "==",scalar(keys %$genes),"\n" if($debug);;
	$cluster_attrs->{'CS3'}++;
      }
      
    }
    else{
	if($alignedstartcount == scalar(keys %$genes)){
	    print "#Class CS2. Inconsistent starts\n" if($debug);;
	    $cluster_attrs->{'CS2'}++;
	}
	else{
	    print "#Class CS3. Unaligned starts ",$alignedstartcount," == ",scalar(keys %$genes),"\n" if($debug);;
	    $cluster_attrs->{'CS3'}++;
	}
    }
    if(scalar(keys %$stops)==1){
      my @stop = keys %$stops;
      if($stops->{$stop[0]}==scalar(keys %$genes)){
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
	    print "#Class CE2. Inconsistent stops\n" if($debug);;
	    $cluster_attrs->{'CE2'}++;
	}
	else{
	    print "#Class CE3. Unaligned stops\n" if($debug);;
	    $cluster_attrs->{'CE3'}++;
	}
	
	#Consider case where we can shift frames and stay open
	
    }
    #Save frequency of starts, stops
    foreach my $feat_name (keys %$genes){
	if(exists $featstarts->{$feat_name}){
	    $feat_attrs->{$feat_name}->{'startfreq='.$starts->{$featstarts->{$feat_name}}}++;
	}
	if(exists $featstops->{$feat_name}){
	    $feat_attrs->{$feat_name}->{'stopfreq='.$stops->{$featstops->{$feat_name}}}++;
	}
    }
    return ($feat_attrs,$cluster_attrs,{'starts'=>$seqstarts,'stops'=>$seqstops});
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
	    if(! exists $codons->{'starts'}->{$seqname}->{$codon}){ #$codon is not annotated on $seqname
		print "#CODON $codon not annotated on $seqname\n" if($debug);;
		#check if $codon is aligned
		#$codon is a tuple of alignment,aligned_column
		my($col,$aln) = split(/\?/,$codon);
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
			    print "#Looking for ORF $start,$end,$orient\n" if($debug);
			    my $seqobj = $db->get_Seq_by_id($seqname);
			    if($seqobj){
				die "Can't find sequence $seqname obj:$seqobj" if(!defined $seqobj);
				my ($neworf,$callorient) = &callORF($seqobj,$start,$end,$orient);
				if(length($neworf)>$MINORF){
				    print "#Calling ORF on strand $callorient start coord = $start\n" if($debug);;
				    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'freq'} = $starts->{$codon};
				    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'neworf'} = $neworf;
				    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'orient'} = $callorient;
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
				    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'start'} = $fmin;
				    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'end'} = $fmax;
				    my($strc,$stpc) = &findCodons($atree,
								  $seqname,
								  $fmin,
								  $fmax,
								  $callorient);
				    if($callorient eq '-'){
					($strc,$stpc) = ($stpc,$strc);
				    }
				    if(ref $strc){
					my($mcol,$align_name) = (@$strc);
					my $startcodon = $mcol.'?'.$align_name;
					#die "Can't find start $mcol,$align_name $callorient,$orient from $seqname $codon" if(!exists $starts->{$startcodon});
					if(!exists $starts->{$startcodon}){
					    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'startfreq'} = 0;
					}
					else{
					    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'startfreq'} = $starts->{$startcodon};
					}
					$codons->{'alt_starts'}->{$seqname}->{$codon}->{'startcol'} = $mcol;
				    }
				    if(ref $stpc){
					my($mcol,$align_name) = (@$stpc);
					my $stopcodon = $mcol.'?'.$align_name;
					#die "Can't find stop $mcol,$align_name $callorient,$orient from $seqname $codon" if(!exists $stops->{$stopcodon});
					if(!exists $stops->{$stopcodon}){
					    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'stopfreq'} = 0;
					}
					else{
					    $codons->{'alt_starts'}->{$seqname}->{$codon}->{'stopfreq'} = $stops->{$stopcodon};					
					}
					$codons->{'alt_starts'}->{$seqname}->{$codon}->{'stopcol'} = $mcol;
				    }
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
    #Report alternative start codons in decreasing order of use
    foreach my $seqname (keys %{$codons->{'alt_starts'}}){
	if(exists $codons->{'alt_starts'}->{$seqname}){
	    foreach my $codon (
			       sort{ #frequency of annotation in the alignment
				   $codons->{'alt_starts'}->{$seqname}->{$b}->{'freq'} <=> $codons->{'alt_starts'}->{$seqname}->{$b}->{$a}->{'freq'}
			       } 
			       keys %{$codons->{'alt_starts'}->{$seqname}}){
		my $althash = $codons->{'alt_starts'}->{$seqname}->{$codon};
		print "#Alt start $seqname start=$codon " if($debug);;
		print "freq=$althash->{'freq'} orient=$althash->{'orient'} start=$althash->{'start'} end=$althash->{'end'}\n" if($debug);;
		if(!exists $seq_attrs->{$seqname}){
		    $seq_attrs->{$seqname} = [];
		}
		if($freq_altstops){ 
		    #check if stop is consistent with at least one other stop
		    if($althash->{'stopfreq'}>1){
			push @{$seq_attrs->{$seqname}},"alt_start=$althash->{'start'}-$althash->{'end'},orient:$althash->{'orient'},len:".($althash->{'end'}-$althash->{'start'}).",startfreq:$althash->{'startfreq'},stopfreq:$althash->{'stopfreq'}";#,startcol:$althash->{'startcol'},stopcol:$althash->{'stopcol'}";
		    }
		}
		else{
		    push @{$seq_attrs->{$seqname}},"alt_start=$althash->{'start'}-$althash->{'end'},orient:$althash->{'orient'},len:".($althash->{'end'}-$althash->{'start'}).",startfreq:$althash->{'startfreq'},stopfreq:$althash->{'stopfreq'}";#,startcol:$althash->{'startcol'},stopcol:$althash->{'stopcol'}";
		}
	    }
	}
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
	    my($startcol,$aln) = split(/\?/,$codon);
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
			    my @indels = &findIndels($atree,$seq,$startcol,$stopcol,$aln,$db);
			    print "#$codon freq:$codons->{'starts'}->{$seq}->{$codon} orient:$orient ORF:$orfstart-$orfend len:",$orfend-$orfstart," aalen:",length($origorf),"\n" if($debug);;
			    foreach my $indel (@indels){
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
				    $neworflist->{$seq}->{$neworfstart.'?'.$neworfend}->{'orf'} = $neworf;
				    $neworflist->{$seq}->{$neworfstart.'?'.$neworfend}->{'start'} = $neworfstart;
				    $neworflist->{$seq}->{$neworfstart.'?'.$neworfend}->{'end'} = $neworfend;
				    $neworflist->{$seq}->{$neworfstart.'?'.$neworfend}->{'orient'} = $neworforient;
				    $neworflist->{$seq}->{$neworfstart.'?'.$neworfend}->{'fsstart'} = $fsstart;
				    $neworflist->{$seq}->{$neworfstart.'?'.$neworfend}->{'orig'} = $origbase;
				    $neworflist->{$seq}->{$neworfstart.'?'.$neworfend}->{'edit'} = $indelbase;
				    $neworflist->{$seq}->{$neworfstart.'?'.$neworfend}->{'startcol'} = $startcol;
				    $neworflist->{$seq}->{$neworfstart.'?'.$neworfend}->{'stopcol'} = $stopcol;
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
	my($col,$aln) = split(/\?/,$codon);
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
    &checkStarts($db,$codons,[keys %$noorfseqs],$seq_attrs);
    return $seq_attrs;
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

sub callORF{
    my($seqobj,$start,$end,$orient) = @_;
    die "Bad start codon $seqobj:$start-$end $orient" if($end < $start);
    my $codon_table = Bio::Tools::CodonTable->new(-id=>11);
    if($seqobj){
	if($orient eq '+'){
            my $seqlen = ($seqobj->length()>$MAXORFLEN) ? $start+$MAXORFLEN : $seqobj->length(); 
	    my $newobj = $seqobj->trunc($start+1,$seqlen);
	    #Check if valid start codon
	    if($codon_table->is_start_codon($newobj->subseq(1,3))){
		my $protein_seq_obj = $newobj->translate(-orf => 1,
							 -codontable_id =>11);
		return ($protein_seq_obj->seq(),$orient);
	    }
	    else{
		print "#callORF trying '-' $seqobj,$start,$end,$orient Bad start codon ",$newobj->subseq(1,3) if($debug);;
		my $seqlen = ($end>$MAXORFLEN) ? $end-$MAXORFLEN : 1;
		my $newobj = $seqobj->trunc($seqlen,$end);
		$newobj = $newobj->revcom();
		#print " REV:",$codon_table->is_start_codon($newobj->subseq(1,3))," ",$newobj->subseq(1,3),"\n";
		if($codon_table->is_start_codon($newobj->subseq(1,3))){
		    my $protein_seq_obj = $newobj->translate(-orf => 1,
							     -codontable_id =>11);
		    
		    return ($protein_seq_obj->seq(),$orient);
		}
		else{
		    print "#WARNING: Skipping callORF $seqobj,$start,$end,$orient. '",$newobj->subseq(1,3),"' is not a valid start codon\n" if($debug);
		}
	    }		
	}
	else{
	    die if($orient ne '-');
            my $seqlen = ($end>$MAXORFLEN) ? $end-$MAXORFLEN : 1;
	    my $newobj = $seqobj->trunc($seqlen,$end);
	    $newobj = $newobj->revcom();
	    #Check if valid start codon
	    if($codon_table->is_start_codon($newobj->subseq(1,3))){
		my $protein_seq_obj = $newobj->translate(-orf => 1,
							 -codontable_id =>11);
		
		return ($protein_seq_obj->seq(),$orient);
	    }
	    else{
		print "#callORF trying '+' $seqobj,$start,$end,$orient Bad start codon ",$newobj->subseq(1,3) if($debug);;
		my $seqlen = ($seqobj->length()>$MAXORFLEN) ? $start+$MAXORFLEN : $seqobj->length(); 
		my $newobj = $seqobj->trunc($start+1,$seqlen);
		if($codon_table->is_start_codon($newobj->subseq(1,3))){
		    my $protein_seq_obj = $newobj->translate(-orf => 1,
							     -codontable_id =>11);
		    
		    return ($protein_seq_obj->seq(),$orient);
		}
		else{
		    print "WARNING: Skipping callORF $seqobj,$start,$end,$orient. '",$newobj->subseq(1,3),"' is not a valid start codon\n";
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
	    print ">CLUSTER_$cluster_id num_seqs=",scalar(keys %$mappedorgs)," num_genes=",scalar(keys %$mappedgenes)," num_alignments=",scalar(@{$mappedgenes->{$query}->{'alignments'}})," classes=$classesstr query=$query alignments=",join(',',@{$mappedgenes->{$query}->{'alignments'}}),"\n";

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
		foreach my $gene (@ogenes){
		    if(exists $feat_attrs->{$gene}){
		      foreach my $c (sort {$a cmp $b} keys %{$feat_attrs->{$gene}}){
		         $classes->{$c}++;
		      }
                    }
		    $longestorf = ($features->{$gene}->[3] > $longestorf) ? $features->{$gene}->[3] : $longestorf;
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
			    #my($orient) = ($alt =~ /orient:(\w+)/);
			    if(!$longer_altstarts || $len>$longestorf){
				push @attrs,$alt;
				push @mappedfeats,[[[$organism,$astart,$aend,'+',($aend-$astart).'M']],"ALTgene$organism$orfidx",'gene'];
				$orfidx++;
				if($organism eq $qseqname){
				    $qfmin = ($astart < $qfmin) ? $astart : $qfmin;	
				    $qfmax = ($aend > $qfmax) ? $aend : $qfmax;
				}
			    }
			    else{
				print "#Skipping $alt $len<$longestorf\n" if($debug);;
			    }
			}
			else{
			    print "#Unknown $alt" if($debug);;
			}
		    }
		}

		foreach my $gene (@ogenes){
		    push @attrs,"feat_orient=$features->{$gene}->[4]";
		    push @attrs,"aln_orient=$mappedgenes->{$gene}->{'relorient'}";
		    if(defined $features->{$gene}->[11]){
			push @attrs,"product=$features->{$gene}->[11]";
		    }
		}
		print join(',',@ogenes),
		"\tWGA$cluster_id",
		"\t$organism",
		"\tcov=",join(',',@ocovs),
		"\tpid=",join(',',@oids),
		"\tqcov=",sprintf("%.2f",$mappedorgs->{$organism}->{'qcov'}/($qfmax-$qfmin)),
		"\t$start-$end",
		"\t",$end-$start,
		"\t",join(';',@attrs),
		"\n";
	    }
	    ##
	    #Report ORFs that are aligned but not annotated
	    foreach my $organism (keys %$new_orfs){
		my $orfidx=0;
		foreach my $alt (@{$new_orfs->{$organism}}){
		    die if(exists $mappedorgs->{$organism});
		    my($astart,$aend) = ($alt =~ /alt_start=(\d+)-(\d+)/);
		    my($len) = ($alt =~ /len:(\d+)/);
		    die "Mismatching lengths $len != $aend - $astart" if($len != ($aend-$astart));
		    #Check that this gene is longer than annotated genes on $organism
		    my @unmappedlist;
		    foreach my $feat_name (keys %$unmappedgenes){
			if($feat2organism->{$feat_name} eq $organism){
			    push @unmappedlist,[$feat_name,$features->{$feat_name}->[3]];
			}
		    }
		    my @longestunmapped = sort {$b->[1] <=> $a->[1]} @unmappedlist;
		    if($len > $longestunmapped[0]){
			push @mappedfeats,[[[$organism,$astart,$aend,'+',($aend-$astart).'M']],"NEWORF$organism$orfidx",'gene'];
			print "NEWORF$organism$orfidx",
			"\tWGA$cluster_id",
			"\t$organism",
			"\tcov=",
			"\tpid=",
			"\t$astart-$aend",
			"\t",$aend-$astart,
			"\t",join(';',@{$new_orfs->{$organism}}),
			"\n";
			$orfidx++;
		    }
		}
	     }
	    if($printalignments){
		print "#Printing query $query $qseqname,$qfmin,$qfmax\n" if($debug);;
		my @isect = $atree->map($qseqname,$qfmin,$qfmax,"alignment");
		#Print all features overlapping the alignment window.
		#This may include addl features than those in the cluster
		my $printedfeats = {};
		foreach my $feat (@isect){
		    my $feat_name = $feat->[0];
		    $feat_name =~ s/gene\://;
		    die "Can't find feature $feat_name" if(!exists $features->{$feat_name});
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
			print "#WARNING Expected gene $feat_name in unmapped list: ".join(',',keys %$unmappedgenes)."\n" if(!exists $unmappedgenes->{$feat_name});
			my $cov = sprintf("c%.1f,i%.1f ",$unmappedgenes->{$feat_name}->{'cov'}/$features->{$feat_name}->[3],
					  $unmappedgenes->{$feat_name}->{'pid'}/$features->{$feat_name}->[3]);
			push @mappedfeats,[[[$seqname,$fmin,$fmax,$orient,($fmax-$fmin).'M']]," $cov *gene:".$feat_name,'gene'];
		    }
		}
		foreach my $align_name (@{$mappedgenes->{$query}->{'alignments'}}){
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
			print "ALIGNMENT:$align_name\n";
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
		    print "#SINGLETON $feat_name\tbest_cluster:$unmapped->{$feat_name}->{'WGA_cluster'}\tcov:";
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
			print "#SINGLETON $feat_name ",join(' ',@$classes);
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
    my $aln_orient=undef;
    if($orient eq '+'){
	if(!$codon_table->is_start_codon($seqobj->subseq($fmin+1,$fmin+2+1))){ #bioperl is 1-base coordinates
	    print "#Bad start codon $fname,$seqname,$fmin,$fmax,$orient codon $fmin+1,$fmin+2+1 ",$seqobj->subseq($fmin+1,$fmin+2+1)," aln_orient:$aln_orient\n" if($debug);
	    return -1;
	} 
	else{
	    #Find start codon + strand
	    $startcodon = &getAlignedCols($atree,$seqname,$fmin,$fmin+3);
	}
	
	if(!$codon_table->is_ter_codon($seqobj->subseq($fmax-3+1,$fmax))){
	    print "#Bad stop $fname,$seqname,$fmin,$fmax,$orient codon $fmax-3+1,$fmax ",$seqobj->subseq($fmax-3+1,$fmax)," aln_orient:$aln_orient\n" if($debug);
	    return -1;
	}
	else{
	    #Find stop codon - strand
	    $stopcodon = &getAlignedCols($atree,$seqname,$fmax-3,$fmax);
	}
    }
    else{
	die "Bad orient $orient" if($orient ne '-');
	if(!$codon_table->is_start_codon(revcom($seqobj->subseq($fmax-3+1,$fmax))->seq())){
	    print "#Bad start codon $fname,$seqname,$fmin,$fmax,$orient codon $fmax-3+1,$fmax ",revcom($seqobj->subseq($fmax-3+1,$fmax))->seq()," aln_orient:$aln_orient\n" if($debug);
	    return -1;
	} 
	else{
	    #Find start codon on - strand
	    $startcodon = &getAlignedCols($atree,$seqname,$fmax-3,$fmax);
	}
	if(!$codon_table->is_ter_codon(revcom($seqobj->subseq($fmin+1,$fmin+3))->seq())){
	    print "#Bad stop codon $fname,$seqname,$fmin,$fmax,$orient codon $fmin+1,$fmin+3 ",revcom($seqobj->subseq($fmin+1,$fmin+3))->seq()," aln_orient:$aln_orient\n" if($debug);
	    return -1;
	} 
	else{
	    #Find stop codon on - strand
	    $stopcodon = &getAlignedCols($atree,$seqname,$fmin,$fmin+3);
	}
    }
    return ($startcodon,$stopcodon);
}


sub getAlignment{
    my($atree,$align_name,$seqname) = @_;
    my $alignment = $atree->{_alignments}->{$align_name}->[0];
    foreach my $i (@$alignment){
	if($i->[0] eq $seqname){
	    return $i;
	}
    }
    print "#Can't find $seqname on alignment $align_name\n" if($verbose);
    return undef;
}

#Look for indels in alignment columns [$codon-$offset,$codon+2]
sub findIndels{
    my($atree,$seq,$startcol,$endcol,$aln,$db) = @_;
    die if($endcol<$startcol);
    print "#Analyzing codon position $startcol in alignment $aln seq $seq \n" if($debug);

    print "#Retrieving alignment matrix for $startcol-$endcol for alignment $aln \n" if($debug);
    my ($mmatrix,$seqmatrix,$names) = $atree->getAlignmentMatrix($aln,$startcol,$endcol,$db);
    print "#Expecting width ",($endcol-$startcol+1)," row count ",scalar(@$mmatrix)," ",scalar(@$names),"\n" if($debug);
    
    my $results = {};
    my @edits;

    my $qryidx;
    my $qrychar;
    my $width;
    for(my $i=0;$i<@$mmatrix;$i++){
	if($names->[$i] eq $seq){
	    $qryidx = $i;
	}
    }
    #Matrix cols start at 0
    for(my $j=0;$j<($endcol-$startcol+1);$j++){
	for(my $i=0;$i<@$mmatrix;$i++){
	    if(substr($mmatrix->[$i],$j,1) ne '.' &&
	       substr($mmatrix->[$i],$j,3) =~ /\./ ){ #gap < length 3
		#column $i has multiple characters
		print "#MUT $i $j ",substr($mmatrix->[$i],$j,1)," $names->[$i] $seq\n" if($debug);
		$results->{$j}++;
	    }

	}
    }
    foreach my $r (keys %$results){
	my $reloffset = $startcol+$r;
	my $freqchar = {};
	for(my $i=0;$i<@$mmatrix;$i++){
	    my $char;
	    if(substr($mmatrix->[$i],$r,1) eq '-'){
		#gap
		$char = substr($mmatrix->[$i],$r,1);
		#die "Unexpected char $i $r $seqmatrix->[$i]->[$r] $mmatrix->[$i]->[$r]" if(defined $seqmatrix->[$i]->[$r]);
	    }
	    else{
		#retrieve base
		#TODO this is slow, improve perf
		$char = substr($seqmatrix->[$i],$r,1);
	    }
	    if($i == $qryidx){
		$qrychar = $char;
	    }
	    die "Bad char '$char'" if(length($char)!=1);
	    $freqchar->{$char}++;
	}
	#report most frequent character
	my @sortedchars = sort {$b <=> $a} (keys %$freqchar);
        #retrieve coordinate on $seq for reloffset
	my $alni = &getAlignment($atree,$aln,$seq);
	my($fsstart,$fsend) = AlignmentTree::columntocoords($alni,$reloffset,$reloffset);
	foreach my $base (@sortedchars){
	    if($base ne $qrychar 
	       #&& $freqchar->{$base}>=$freqchar->{$qrychar}	    #only consider bases that occur more frequently than 
	       #&& $freqchar->{$base}>=scalar(@$mmatrix)/2){  	    #optionally also in majority of sequences
	       ){
		print "#ALT col:$reloffset coord:$fsstart-$fsend base:$base freq:$freqchar->{$base} $seq:$qrychar $freqchar->{$qrychar}\n" if($debug);
		push @edits,[$fsstart,$base,$qrychar,$reloffset];
	    }
	    else{
		#last;#can shortcircuit, only consider more frequent bases
	    }
	}
    }
    return @edits;
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
