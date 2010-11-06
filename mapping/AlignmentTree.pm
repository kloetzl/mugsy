package AlignmentTree;

#AlignedIntervalTree is an interval tree with the additions that
#stored intervals 1) may contain a correspondence map, such as an
#alignment and 2) can be oriented for DNA sequences.

#The data structure supports retrieval of corresponding,aligned intervals
#The data structure also supports discontinuous intervals

#Each interval in the structure is associated with a single coordinate
#system or sequence and has an orientation '+','-'

#The data structure used to represent an interval and an alignment is
#[[seqname1,start1,end1,orientation2,cigarstring1,tag1_0,...tag1_N],
# [seqname2,start2,end2,orientation2,cigarstring2,tag2_0,...tag2_N],...,]

#Represented in the code as $alignobj = [$alni_1,$alni_2];

#
#insert(interval/alignment)     - insert an interval or alignment (a series of mapped intervals)
#find(seq,start,end)            - retrieve intervals that overlap start,end on seq
#intersect(seq,start,end)       - retrieve corresponding,aligned intervals that overlap start,end on seq
#                                 
#map(seq,start,end)             - retrieve intervals that overlap the range specified by any intersecting intervals 
#                                 intersecting intervals are obtained if there exists an alignment that spans start,end on seq
# Definitions
# interval
# alignment -  a series of mapped intervals

use strict;
use Math::Random qw(random_uniform);
use POSIX qw(ceil floor);
use IntervalTree;
use Bit::Vector;
use Storable qw(store retrieve);

#remove only using for translation machinery and revcom
use Bio::Perl;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::Tools::CodonTable;

$Storable::Deparse = 1;
$Storable::Eval = 1;
my $DEBUG=0;
my $BITV_SIZE=10000000; #10MB largest single aligned region

#my $aligntoken="WGA";
my $aligntoken="WGA";

sub new{
    my $classname = shift;
    my $self = {};
    bless($self,$classname);
    $self->{_itrees} = {};
    $self->{_alignments} = {}; #Saved as [alignref,bitvector,align_width]
    #Support for filtering output using a phylogenetic profile of genomes
    #Implemented using bitmasks
    $self->{_maxbits} = 1000;

    $self->{_doremoveoverlaps}=0;
    $self->{_bits} = 0;
    $self->{_bitlookup} = {};
    $self->{_bitmask} = new Bit::Vector($self->{_maxbits});
    $self->{_defaultmask} = new Bit::Vector($self->{_maxbits});
    $self->{_debug}=$DEBUG;
    return $self;
}

sub serialize{
    my($self,$file) = @_;
    $self->{_bitmaskstr} = $self->{_bitmask}->to_Enum();
    $self->{_defaultmaskstr} = $self->{_bitmask}->to_Enum();
    return Storable::store($self,$file);
}
sub deserialize{
    my($file) = @_;
    my $atree = Storable::retrieve($file);
    $atree->{_bitmask} = new Bit::Vector($atree->{_maxbits});
    $atree->{_bitmask}->from_Enum($atree->{_bitmaskstr});
    $atree->{_defaultmask} = new Bit::Vector($atree->{_maxbits});
    $atree->{_defaultmask}->from_Enum($atree->{_bitmaskstr});
    return $atree;
}

#Require output contains one or more tags
#An example of a tag is a genome name
sub filter{
    my($self) = shift;
    foreach my $tag (@_){
	if(!exists $self->{_bitlookup}->{$tag}){
	    $self->{_bitlookup}->{$tag} = $self->{_bits}++;
	}
	$self->{_bitmask}->flip($self->{_bitlookup}->{$tag});
    }
    $self->{_bitmask}->Union($self->{_defaultmask},$self->{_bitmask});
}

sub clear_filter{
    my($self) = shift;
    $self->{_bitmask}->Empty();
}

#
#Insert an interval or alignment
#insert(
#       {[seqname,start,end,orientation,cigarstring,tag0,...tagN]},
#       uniquename,
#       tags
#      )
#A unique identifier for the alignment, uniquename, must be provided
#seqname - must be a uniquename for the coordinate system containing interval [start,end]
#start - beginning of interval 0-based
#end - end of interval 0-based
#orientation - 2
#cigarstring - in the UCSC format (#M#S#I#D#X) indicates the continuity of the alignment over the interval
#tag0...tagN - zero or more tags that can be used by filtering functions. Tags can be specified on either the alignment or the interval

#Intervals and alignments are stored in a consistent manner. An
#alignment is a set of intervals with a correspondence map. Single
#annotated intervals, like genes, are stored as an alignment alignment
#with only a single interval.  The correspondence map is an identity
#map in this case.

sub insert{
    my($self,$alignmentref,$name,@tags) = @_;
    my $genomelookup = {};
    my $alignment_bv = new Bit::Vector(1000);
    my $align_width = 0;
    die "Bad alignment passed to insert($alignmentref). Alignment needs to be a ref to an array" if(!ref($alignmentref));
    foreach my $align (@$alignmentref){
	die "Bad alignment passed to insert($align). Alignment needs to be a ref to an array" if(!ref($align));
	#print "INSERTING ",join(',',@$align),"\n";
	my $seqname = $align->[0];
	my $start = $align->[1];
	my $end = $align->[2];
	my $orientation = $align->[3];
	if($align->[4]){
	    #Check that column count is consistent
	    my ($cigs,$columncount) = &get_cigs($align->[4]);
	    $align_width = $columncount if(!$align_width);
	    if($columncount != $align_width){
		&printAlignmentDebug($alignmentref,\*STDERR);
		die "Bad input. Mismatched column count $columncount in $align->[4], expecting $align_width";
	    }
	}
	if($orientation =~ /\d/){
	    if($orientation>0){
		$orientation = '+';
	    }
	    else{
		$orientation = '-';
	    }
	} 
	$align->[3]=$orientation;
	die "Bad orient $orientation ".join(',',@$align)."\n" if($orientation ne '-' && $orientation ne '+');
	#Store tags in bit vector
	for(my $i=5;$i<@$align;$i++){
	    my $tag = $align->[$i];
	    if(!exists $self->{_bitlookup}->{$tag}){
		$self->{_bitlookup}->{$tag} = $self->{_bits}++;
	    }
	    $alignment_bv->Bit_On($self->{_bitlookup}->{$tag});
	}
	if(!exists $self->{_itrees}->{$seqname}){
	    $self->{_itrees}->{$seqname} = new IntervalTree($start,$end,$name,$orientation);
	}
	else{
	    $self->{_itrees}->{$seqname}->insert($start,$end,$name,$orientation);
	}
    }

    #Store tags in bit vector
    foreach my $tag (@tags){
	if(!exists $self->{_bitlookup}->{$tag}){
	    $self->{_bitlookup}->{$tag} = $self->{_bits}++;
	}
	#print STDERR "Adding tag $tag on $self->{_bitlookup}->{$tag} $self->{_defaultmask}\n";
	$alignment_bv->Bit_On($self->{_bitlookup}->{$tag});
	$self->{_defaultmask}->Bit_On($self->{_bitlookup}->{$tag});
    }
#    print "Masks ",$self->{_defaultmask}->Norm()," ",$self->{_bitmask}->Norm(),"\n";
    if(0 && exists $self->{_alignments}->{$name}){
	print STDERR "Duplicate feature $name already stored. Skipping this one\n";
    }
    else{
	$self->{_alignments}->{$name} = [$alignmentref,$alignment_bv,$align_width];
    }
}

#
#Find all intersecting alignments in interval (query.start,query.end) from query.seqname 
#Returns (start,end) coordinates on seqname of all matching alignments
#
#intersect(query.seqname,query.start,query.end,tags)

#returns [alignname,seqname,start,end,coverage,pid,queryorient,matchorient]
#0-alignname
#1-seqname
#2-start
#3-end
#4-coverage is number of corresponding characters between start,end
#5-pid is number of identical characters between start,end
#6-queryorient is orientation of the matching aligned query interval query.seqname:query.start-query.end 
#7-matchorient is orientation of the matching aligned interval seqname:start-end

sub intersect{
    my($self,$qseqname,$qstart,$qend,@qtags) = @_;
    my @results;
    #$self->filter(@qtags);
    if(exists $self->{_itrees}->{$qseqname}){
	print "Querying $qseqname:$qstart,$qend with qtags $qtags[0]\n" if($self->{_debug});
	#(1) Find all intersecting features on [$qstart,$qend]
	#returns IntervalTree::intersect returns an array of interval names
	my @alignments = $self->{_itrees}->{$qseqname}->intersect($qstart,$qend); 
	#Optionally remove fully nested intervals
	if($self->{_doremoveoverlaps}){
	    @alignments = $self->removeOverlaps(\@alignments,$qseqname);
	}
	foreach my $align_name (@alignments){
	    die "Overlapping interval $align_name not found" if(! exists $self->{_alignments}->{$align_name});
	    my($alignobj,$alignment_bv,$align_width) = @{$self->{_alignments}->{$align_name}};
	    if($align_name =~ /$qtags[0]/){
		print "Overlapping feature $align_name\n" if($self->{_debug});
		print "MATCH $align_name query:$qseqname $qstart-$qend . Number of seqs ",scalar(@$alignobj),"\n" if($self->{_debug});
		#(2) Crop interval [$qstart,$qend] to the alignment 
		my ($qmstart,$qmend,$queryorient) = &matchinginterval($alignobj,$qseqname,$qstart,$qend);
		if(!defined $qmstart || !defined $qmend){
		    #Error condition
		    print "WARNING. print unexpected overlapping alignments for query $qstart,$qend\n";
		    foreach my $align_name2 (@alignments){
			my($alignobj2,$alignment_bv2,$align_width2) = @{$self->{_alignments}->{$align_name2}};
			foreach my $alni2 (@$alignobj2){
			    print "$align_name2  $alignobj2 ",join(' ',@$alni2),"\n";
			}
		    }
		    die "Bad overlapping alignments";
		}
		die if($qmstart<$qstart);
		die if($qmend>$qend);
		die "Invalid matching interval coords:$qmend-$qmstart from query $qstart-$qend\n" if($qmend<=$qmstart);
		my $queryspancheck=0;
		if($qstart == $qmstart && $qend == $qmend){
		    print "Alignment fully spans query\n" if($self->{_debug});
		}
		else{
		    if($qstart != $qmstart && $qend != $qmend){
			print "Query fully spans alignment\n" if($self->{_debug});
			$queryspancheck=1;
		    }
		}
		print "ISECT: $align_name QUERY:$qseqname $qstart-$qend mapped:$qmstart-$qmend len:",$qmend-$qmstart,"\n" if($self->{_debug});
		#(3) Convert from genomic coords to alignment column. 1->alignment_width
		my ($qcolumnstart,$qcolumnend,$querybv) = &coordstocolumn($alignobj,$qseqname,$qmstart,$qmend);
		#$querybv stores a bitmatrix from $qcolumnstart-$qcolumnend indicating if sequence $seqname is aligned in the interval
		print "MAPPED $qseqname:$qmstart-$qmend len:",$qmend-$qmstart,
		" to column coords $qseqname:$qcolumnstart-$qcolumnend len:",$qcolumnend-$qcolumnstart+1,"\n" if($self->{_debug});
		print "Transform $qseqname:$qmstart-$qmend to column coords $qcolumnstart-$qcolumnend\n" if($self->{_debug});
		
		foreach my $alni (@$alignobj){
		    die "Invalid zero length matching interval $alni->[2]-$alni->[1]\n" if($alni->[2] - $alni->[1]<=0);
		    print "Converting $align_name $alni->[0]:$alni->[1]-$alni->[2] from $qseqname:$qmstart-$qmend using column coords $qcolumnstart-$qcolumnend\n" if($self->{_debug});
		    my ($mseq,$malign_start,$malign_end,$morient) = @$alni;
		    print "ALNI: $mseq,$malign_start,$malign_end,$morient\n" if($self->{_debug});
		    die if(@$alni>5);
		    
		    #(4) Crop aligned feature. Convert back from alignment column to genomic coords on $mseq
		    #$currbv stores a bitmatrix from $qcolumnstart-$qcolumnend indicating if sequence $mseq is aligned in the interval
		    my($s,$e,$currbv) = &columntocoords($alni,$qcolumnstart,$qcolumnend,$querybv);
		    if($mseq eq $qseqname){
			die if($s != $qmstart);
			die if($e != $qmend);
			die if($morient ne $queryorient);
		    }
		    
		    #Check the actual number of aligned columns
		    my $pid=0;
		    if($mseq eq $qseqname){
			my ($qs1,$qe1) = &coordstocolumn($alignobj,$qseqname,$s,$e);
			print "Checking for matching characters between col:$qcolumnstart-$qcolumnend $qs1-$qe1 coords:$s-$e\n" if($self->{_debug});
		    }

		    my $intersectbv = new Bit::Vector($querybv->Size());
		    $intersectbv->Intersection($querybv,$currbv); 
		    for(my $i=$qcolumnstart;$i<=$qcolumnend;$i++){
			if($intersectbv->bit_test($i)==1){
			    $pid++;
			}
		    }
		    print "Intersect matches: ",$intersectbv->Norm(),"\n" if($self->{_debug});
		    print "Query matches: ",$querybv->Norm(),"\n" if($self->{_debug});
		    print "Matches in the interval $qcolumnstart-$qcolumnend:$pid\n" if($self->{_debug});

		    if($e-$s>0){
			die "Bad pid" if($pid==0);
			die "Invalid zero length matching interval $e-$s\n" if($e-$s<=0);
			print "($qcolumnend-$qcolumnstart) - ($qmend-$qmstart)\n" if($self->{_debug});
			my $numgaps_query = ($qcolumnend-$qcolumnstart+1) - ($qmend-$qmstart);
			my $querypid = ($qcolumnend-$qcolumnstart+1) - $numgaps_query;
			my $coverage = $e-$s;
			print "Num_query_gaps=$numgaps_query\nNum_qry_matches=$querypid\nNum_hit_matches=$pid\n" if($self->{_debug});
			die "$pid<1" if($pid<1);
			die "$pid>$coverage" if($pid>$coverage);
			print "RESULT: $align_name,$mseq,$s,$e,$coverage,$pid\n" if($self->{_debug});
			#Intersect result is a $alni,$coverage,$pid
			push @results,[$align_name,$mseq,$s,$e,$coverage,$pid,$queryorient,$morient];
		    }
		    else{
			#Entirely contained within a gap
			die "Bad pid" if($pid !=0);
			die "Bad coordinates $align_name,$mseq,$s,$e\n" if($e<$s);
			print "NORES: Skipping $align_name,$mseq,$s,$e. Mapped in a gap\n" if($self->{_debug});
		    }
		}
	    }
	    else{
		print "Skipping $align_name does not match $qtags[0]\n" if($self->{_debug});;
	    }
	}
    }
    else{
	#Skip, nothing to find
	if($self->{_debug}){
	    print "Interal trees for ",join(',',keys %{$self->{_itrees}}),"\n";
	    print "No interval tree for sequence [$qseqname] $self->{_itrees}->{$qseqname}\n";
	}
    }
    #$self->clear_filter();
    return @results;
}
#
#Map a coordinate (query.start,query.end) from query.seqname to
#intersecting alignments on match_i.seqname..match_j.seqname with
#coordinates (match_i.start,match_i.end...match_j.start,match_j.end)
#
#map(seqname,start,end,type)
#returns [match_name,$mseq,$mstart,$mend,$mcoverage,align_name,seq_name,$coveraged]
sub map{
    my($self,$qseqname,$qstart,$qend,@qtags) = @_;
    my @results;
    if(exists $self->{_itrees}->{$qseqname}){
	print "Finding WGA alignments on $qseqname,$qstart,$qend\n" if($self->{_debug});
	#(1)Retrieve all the alignments on genomic coords qstart-qend
	$self->{_doremoveoverlaps}=1;
	my @isects = $self->intersect($qseqname,$qstart,$qend,$aligntoken);
	$self->{_doremoveoverlaps}=0;
	print "FOUND ", scalar(@isects)," alignments\n" if($self->{_debug});

	#Currently assuming non overlapping alignments and the total
	#coverage,pid must be less than the query length $qend-$qstart

	my $totalqcoverage=0;
	my $totalqid=0;
	my $qcoverage=undef;
	my $qmstart=undef;
	my $qmend=undef;
	my $qmorient=undef;

	#(2)Determine the min-max spanning interval over all matching alignments to the query
	#intersect() already provides the query interval [$qstart,$qend] crop to the overlapping alignment(s)
	foreach my $isectn (@isects){
	    my($align_name,$seq,$start,$end,$coverage,$pid,$qorient,$orient) = @$isectn;
	    print "Looking for $qseqname in $align_name,$seq,$start,$end,$coverage,$pid\n" if($self->{_debug}); 
	    if($seq eq $qseqname){
		die "Mismatched orient $qmorient != $orient" if($qorient ne $orient);
		die "$end>$qend" if($end>$qend);
		die "$start<$qstart" if($start<$qstart);
		if(defined $qmstart || defined $qmend){
		    print "#Duplicate $seq already found in $align_name. Multiple alignments spanning query\n" if($self->{_debug}); 
		    $qcoverage=$coverage;
		    $totalqcoverage+=$coverage;
		    $totalqid+=$pid;
		    $qmstart=$start<$qmstart ? $start : $qmstart;
		    $qmend=$end>$qmend ? $end : $qmend;
		    if(defined $qmorient && $orient ne $qmorient){
			#print "WARNING multiple matching alignments to $qseqname,$qstart,$qend with inconsistent orientations.  $align_name:$orient ne $qmorient\n";
			$qmorient='?';
		    }
		    else{
			$qmorient=$orient;
		    }

		}
		else{
		    die if(defined $qmorient);
		    $qcoverage=$coverage;
		    $totalqcoverage+=$coverage;
		    $totalqid+=$pid;
		    $qmstart=$start;
		    $qmend=$end;
		    $qmorient=$orient;
		}
	    }
	}
	
	#
	#(3)Map the spanning interval [$qmstart,$qmend] to the rest of the sequences in the alignment
	
	foreach my $isectn (@isects){
	    die if(!defined $qcoverage);
	    die if(!defined $qmstart || !defined $qmend);
	    my($align_name,$seq,$start,$end,$coverage,$pid,$qaln_orient,$aln_orient) = @$isectn;
	    my($alignobj,$alignment_bv,$align_width) = $self->getAlignment($align_name);
	    my($qfstart,$qfend,$qforient) = &matchinginterval($alignobj,$qseqname,$qmstart,$qmend);
	    my ($qfscolumnstart,$qfscolumnend,$fsquerybv) = &coordstocolumn($alignobj,$qseqname,$qfstart,$qfend,$qforient);
	    #print "#Mapping with alignment $align_name $seq $start-$end cov:$coverage,pid:$pid,qaln_orient:$qaln_orient,aln_orient:$aln_orient\n";
	    #Query coverage should correspond to interval start-end
	    die if($coverage != ($end-$start));
	    die "$qaln_orient ne $qmorient" if($qmorient ne '?' && $qaln_orient ne $qmorient);
	    print "$align_name:$seq $start,$end ",$end-$start," query_start:$qstart query_end:$qend query_coverage:$qcoverage ",$qcoverage/($qend-$qstart),"\n" if($self->{_debug});
	    #(4)Find features in the mapped interval
	    my @misects = $self->intersect($seq,$start,$end,"gene");
	    foreach my $fisectn (@misects){
		my($fname,$fseq,$fstart,$fend,$fcoverage,$fpid,$forient1,$forient2) = @$fisectn;

		#Need intersection of $start,$fend $qmstart,$qmend to get proper $pid and $cov
		my ($fscolumnstart,$fscolumnend,$fsbv) = &coordstocolumn($alignobj,$seq,$fstart,$fend);
		#die "$fscolumnstart != $qfscolumnstart" if($fscolumnstart != $qfscolumnstart);
		#die "$fscolumnend != $qfscolumnend" if($fscolumnend != $qfscolumnend);
		my $ipid=0;
		my $intersectbv = new Bit::Vector($fsquerybv->Size());
		$intersectbv->Intersection($fsquerybv,$fsbv); 
		for(my $i=$fscolumnstart;$i<=$fscolumnend;$i++){
		    if($intersectbv->bit_test($i)==1){
			$ipid++;
		    }
		}
		die "Bad number of matching columns $ipid>($fend-$fstart) $fscolumnstart-$fscolumnend" if($ipid>($fend-$fstart));
		die if($forient1 ne $forient2);
		die if($seq ne $fseq);
		if($fseq eq $qseqname){
		    die "$fstart<$qmstart query:$seq,$start,$end $fname,$fseq,$fstart,$fend,$fcoverage,$fpid" if($fstart<$qmstart);
		    die "$fend>$qmend query:$seq,$start,$end $fname,$fseq,$fstart,$fend,$fcoverage,$fpid" if($fend>$qmend);
		}
		print "Adding result $fname,$fseq,$fstart,$fend,$fcoverage,$align_name,$seq,$qcoverage,$fpid\n" if($self->{_debug});
		#push @results,[$fname,$fseq,$fstart,$fend,$coverage,$align_name,$seq,$fcoverage,$ipid,$isectn,$qaln_orient,$aln_orient,$forient1];
		#Determine span on query
		my ($qfsstart,$qfsend) = &columntocoords($self->getAlignedInterval($align_name,$qseqname),$fscolumnstart,$fscolumnend);
		push @results,[$fname,$fseq,$fstart,$fend,$qfsend-$qfsstart,$align_name,$seq,$fcoverage,$ipid,$isectn,$qaln_orient,$aln_orient,$forient1];
	    }
	    print "Finished mapping alignment $align_name\n" if($self->{_debug});
	}
	if($totalqcoverage>($qend-$qstart)){
#TODO die on bad coverage die "Bad coverage $totalqcoverage>($qend-$qstart) ".($qend-$qstart) if($totalqcoverage>($qend-$qstart));
	    $totalqcoverage=($qend-$qstart);	  
	}
	if($totalqid>($qend-$qstart)){
#TODO die "Bad identity ($totalqid>($qend-$qstart) " if($totalqid>($qend-$qstart));
	    $totalqid = ($qend-$qstart);
	}
    } 
    return @results;
}

sub matchinginterval{
    my($alignobj,$qseqname,$qstart,$qend) = @_;
    my $start=undef;
    my $end=undef;
    my $orient=undef;
    foreach my $alni (@$alignobj){
	if($alni->[0] eq $qseqname){
	    print "HIT on $alni->[0] query=$qseqname:$qstart-$qend ; interval=$qseqname:$alni->[1]-$alni->[2] $alni->[3]\n" if($DEBUG);
	    if(($qstart < $alni->[1] && $qend < $alni->[1]) || ($qstart > $alni->[2] && $qend > $alni->[2])){
		print "WARNING: Invalid matching interval. Alignment interval $alni->[0]:$alni->[1]-$alni->[2] not contained in interval $qseqname:$qstart-$qend\n";
		&printAlignmentDebug($alignobj);
		return ($start,$end,$orient);
	    }

	    $start = $qstart < $alni->[1] ? $alni->[1] : $qstart;
	    $end = $qend < $alni->[2] ? $qend : $alni->[2];
	    if(defined $orient && $orient ne $alni->[3]){
		print "WARNING multiple matching alignments to $qseqname,$qstart,$qend with inconsistent orientations. $orient ne $alni->[3]\n";
		&printAlignmentDebug($alignobj);
		die "Multiple copies of a sequence per alignment not supported";
		#Keep existing orient
	    }
	    else{
		$orient = $alni->[3];
	    }
	}
	else{
	    print "Checked $alni->[0] in obj size ",scalar(@$alignobj),"\n" if($DEBUG);
	}
    }
    die "Cannot find $qseqname,$qstart,$qend in alignment. Returned ($start,$end,$orient)" if(!defined $start || !defined $end);
    return ($start,$end,$orient);
}

#Determine column indices in an alignment matrix
#that correspond to interval $coord1-$coord2 on sequence $qseqname
#The alignment matrix is assumed to have one sequence per row with no
#sequence appearing more than oncequerybv
#The format is
#@$alignobj = [seqname,start,end,orient,cigar]
#start<end is specified in 0 start interbase coordinates
#cigar specifies the continuity of the alignment of the interval
#when orient=='-' the cigar string specifies the alignment end->start
#otherwise the cigar string species the alignment start->end
#Column coordinates are 1 start, numbering bases/columns in an alignment matrix
#Genomic coordinates are 0 start, interbase
sub coordstocolumn{
    my($alignobj,$qseqname,$coord1,$coord2) = @_;
    die "Expecting 0 start, interbase coordinates $coord1<=$coord2" if($coord1>=$coord2);
    #Column position in the alignment
    #Starting from position 1
    my $columnstart;
    my $columnend;
    #Bit vector keeps track of aligned columns
    #Starting at column 1, column index 0 is ignored
    my $querybv = new Bit::Vector($BITV_SIZE); #setting max aligned interval at 10MB
    foreach my $alni (@$alignobj){
	if($alni->[0] eq $qseqname){

	    if($coord1<$alni->[1] || $coord1>$alni->[2]){
		&printAlignmentDebug($alignobj);
		die "Start position $coord1 is not contained in interval $qseqname:$alni->[1]-$alni->[2]";
	    }
	    if($coord2<$alni->[1] || $coord2>$alni->[2]){
		&printAlignmentDebug($alignobj);
		die "End position $coord2 is not contained in interval $qseqname:$alni->[1]-$alni->[2]";
	    }

	    my $offsetstart;
	    my $offsetend;
	    my $orient = $alni->[3];
	    print "COORDSTOCOLUMN for seq $qseqname and coordinates $coord1-$coord2 and orient $orient\n" if($DEBUG);

	    die "Bad orient: $orient\n" if($orient ne '-' && $orient ne '+');
	    #$alignobj is a collinear segment
	    #Determing offsets into the aligned segment
	    #$alni->[1] |-----------------| $alni->[2]
	    #coord1        |------|      coord2
	    #+ strand
	    #        ------|             offsetstart
	    #        -------------|      offsetend
	    #Cigar string is relative to $alni->[1] --> $alni->[1]
	    #Columnstart,end are relative to $alni->[1]              

	    #- strand
	    #                     |----- offsetstart
	    #              |------------ offsetend
	    #Cigar string is relative to $alni->[2] --> $alni->[1]      
	    #Columnstart,end are relative to $alni->[2]                             
	    if($orient eq '-'){
		#offset is from end ($coord2) of segment match
		# A T G C A T
		#0 1 2 3 4 5 6
		#  |---------|   i->[1] -> i[2] == 1-6 
		#    |-----|     coord1 -> coord2 == 2-5
		#         *      offsetstart = 6-5+1 == 2; column 2 in the alignment matrix
		# X X M M M X    alignment matrix
		#     *          offsetend = 6-2 == 4; column 4 in the alignment matrix
		#0 1 2 3 4 5 6
		# A T G C A T
		$offsetstart = $alni->[2] - $coord2+1;
		$offsetend = $alni->[2] - $coord1;
		die "$offsetstart = $alni->[2] - $coord2+1; $offsetend = $alni->[2] - $coord1;" if($offsetend < $offsetstart);
	    }
	    else{
		#offset is from beginning ($coord1) of segment match
		$offsetstart = $coord1 - $alni->[1]+1;
		$offsetend = $coord2 - $alni->[1];
		die "$offsetstart < $offsetend " if($offsetend < $offsetstart);
	    }
	    die if($offsetstart<1);
	    die if($offsetend>($alni->[2]-$alni->[1]));
	    #If cigar string
	    if($alni->[4]){
		print "coordstocolumn using cig $alni->[0] $alni->[4]\n" if($DEBUG);
		print "Looking for offset $offsetstart-$offsetend:$orient in match $alni->[1]-$alni->[2] of length ",
		$alni->[2]-$alni->[1]," orient:$orient\n" if($DEBUG);

		my ($cigs,$columncount) = &get_cigs($alni->[4]);
		die "$offsetend>$columncount. Check cigar string $alni->[4], appears to be incorrect length for $coord1-$coord2" if($offsetend>$columncount);
		my $currcount2=0;
		foreach my $c2 (@$cigs){
		    my($count2,$char2) = @$c2;
		    if($char2 eq 'M'){
			#|*      |$currcount2
			#| *     |$currcount2+1
			#|     * |$currcount2+$count2
			#| MMMXX |3M,$count2==3
			#| 11100 |bitvector
			$querybv->Interval_Fill($currcount2+1,$currcount2+$count2);
		    }
		    $currcount2+=$count2;
		}

		my $matches=0;
		my $currcount=0;

		my $foundstart=0;
		my $foundend=0;			
		
		foreach my $c (@$cigs){
		    my($count,$char) = @$c;
		    if($char eq 'M'){
			#0 1 2 3 4 5 6
			#  |---------|   i->[1] -> i[2] == 1-6 
			#     *          offsetstart
			#         *      offsetend
			# X X M M M X    alignment matrix
			#   *            matches
			#0 1 2 3 4 5 6
			# A T G C A T
			if($count+$matches>=$offsetstart){
			    if(!$foundstart){
				print "FOUNDSTART $currcount $matches $count$char\n" if($DEBUG>1);
				$columnstart = $currcount+($offsetstart-$matches);
				print "columnstart=$columnstart\n" if($DEBUG>1);
				$foundstart=1;
			    }
			}
			if($count+$matches>=$offsetend){
			    if(!$foundend){
				print "FOUNDEND $currcount $matches $count$char\n" if($DEBUG>1);
				#$columnend=$currcount+($offsetend-$matches)-1;
				$columnend=$currcount+($offsetend-$matches);
				print "columnend=$columnend\n" if($DEBUG>1);
				$foundend=1;
				last;
			    }
			}
			$matches+=$count;
			$currcount+=$count;
		    }
		    elsif($char eq 'X'){
			$currcount+=$count;
		    }
		}
		die "Could not find start or end " if(!$foundstart || !$foundend);
		die if($columnstart<1);
		die if($columnend>$columncount);
	    }
	    else{
		#No cigar string 
		#Assume interval aligns at its entire length
		$columnstart=$offsetstart;
		$columnend=$offsetend;
		$querybv->Interval_Fill($columnstart,$columnend);
	    }
	    last;
	}
    }
    die "Can't map $columnstart-$columnend" if(!defined $columnstart || !defined $columnend);
    return ($columnstart,$columnend,$querybv);
}

#Maps the alignment columns $columnstart-$columnend to genomic coordinates
#In the case where the specified columns map to gaps, the coordinate corresponding to the 
#next matching column is returned 
#Column coordinates are 1 start, numbering bases,gaps/columns in an alignment matrix
#Genomic coordinates are 0 start, interbase
sub columntocoords{
    my($aln,$columnstart,$columnend,$querybv) = @_;
    die "Columnstart $columnstart must be >= 1. aln:$aln" if($columnstart<1);
    #columnstart is start relative to beginning of alignment matrix
    #columnend is end relative to beginning of alignment matrix
    #|MMMIIMMMIIIMMMXXXMM|
    # -----|               columnstart
    # ------------|        columnend
    #Alignment $aln is in genomic coordinates
    #start,end are mapped genomic coordinates for aln relative to offsets specified by $columnstart,$columnend
    my $start; 
    my $end;
    #Bit vector keeps track of aligned columns
    #Starting at column 1
    if(! defined $querybv){
	$querybv = Bit::Vector->new($BITV_SIZE);
	$querybv->Interval_Fill($columnstart,$columnend);
    }
    my $currbv = Bit::Vector->new($querybv->Size()); 
    my $orient = $aln->[3];
    die "Bad orient $orient\n" if($orient ne '-' && $orient ne '+');
    #Check if there is a cigar string
    print "Converting for alni with orient $orient\n" if($DEBUG);
    if(length($aln->[4])>0){
	print "columntocoords using cig $aln->[0] $aln->[4]\n" if($DEBUG > 1);
	my ($cigs,$count) = &get_cigs($aln->[4]);
	die "Columnstart $columnstart must be >= 1. aln:$aln" if($columnstart<1);
	die "Columnend $columnend > cigar count $count. Check cigar string $aln->[4] for $aln->[1]-$aln->[2]" if($columnend>$count);
	print "Looking for offset $columnstart-$columnend:$orient in match $aln->[1]-$aln->[2] of length ",$aln->[2]-$aln->[1]," orient:$orient\n" if($DEBUG > 1);
		
	my $matches=0;
	my $currcount=0;
	
	my $foundstart=0;
	my $foundend=0;
	my $startcount;
	my $currcount2=0;
	foreach my $c2 (@$cigs){
	    my($count2,$char2) = @$c2;
	    if($char2 eq 'M'){
		$currbv->Interval_Fill($currcount2+1,$currcount2+$count2);
	    }
	    $currcount2+=$count2;
	}
	foreach my $c (@$cigs){
	    my($count,$char) = @$c;
	    print "Analyzing $count$char\n" if($DEBUG>1);
	    if($count+$currcount>=$columnstart){
		if(!$foundstart){
		    print "FOUNDSTART $currcount $matches $count$char\n" if($DEBUG > 1);
		    if($char eq 'M'){
			print "START IN MATCH $orient\n" if($DEBUG > 1);
			#|----*------ columnstart
			# --|         currcount
			#    MMMM     count        columnstart-currcount=number of contributing matches in current cig
			if($orient eq '-'){
			    $start = $aln->[2]-$matches-($columnstart-$currcount)+1;
			}
			else{
			    #Start is alignment start (s1) + matching columns + number of matches in current cigar
			    $start = $aln->[1]+$matches+($columnstart-$currcount)-1;
			}
		    }
		    else{
                        #Report NEXT matching position
			print "START IN GAP $orient\n" if($DEBUG > 1);
			#No match to $->[0]; in gap
			#Report next matching position
			#in the alignment between query and current sequence
			#(the next matching position past the gap)

                        #|----*------ columnstart
			# --|         currcount
			#    XXXX     count

			#Need to account for case where
			#next matching position in current sequence is a gap in the query
			#Use the bit vectors to find next matching position between current seq and query
			my $intersectbv = new Bit::Vector($querybv->Size());
			$intersectbv->Intersection($querybv,$currbv); 
			die if($currbv->bit_test($currcount+1));
			die if($intersectbv->bit_test($currcount+1));
			my($imin,$imax) = $intersectbv->Interval_Scan_inc($currcount+1);
			my($cmin,$cmax) = $currbv->Interval_Scan_inc($currcount+1);
			if(! defined $cmin || ! defined $imin){
			    print "INGAP Can't find matching position in $aln->[0] > columnstart:$columnstart. Returning no mapping\n" if($DEBUG > 1);
			    return ($aln->[1],$aln->[1],$currbv);
			}
			if($imin>$columnend){
			    #|-------*MMMMMMM*-----------| Query columnstart -> columnend
			    #        *       *             columnstart,columnend
			    #                    *         cmin>columnend
			    #|----*XXXXXXXXXXXXXXMM------| 
			    #  
			    print "INGAP Alignment occurs entirely within a gap $cmin>$columnend\n" if($DEBUG > 1);
			    return ($aln->[1],$aln->[1],$currbv);
			}
			die "$imin,$imax $cmin,$cmax" if($imin<$cmin);
			#num matching bits set between 
			#$imin-$cmin is number of matches in current seq until the 
			#next matching position in the query
			my $nummatches=0;
			if($cmin>0){
			    die "Bad match index after scan $cmin" if($currbv->bit_test($cmin)!=1);
			}
			for(my $i=$cmin;$i<=$imin;$i++){
			    if($currbv->bit_test($i)){
				$nummatches++;
			    }
			}
			die if($cmin==0 && $nummatches!=0);
			print "INGAP $nummatches matches until next match between query and $aln->[0] between columns $cmin-$imin\n" if($DEBUG > 1);
			if($orient eq '-'){
			    $start = $aln->[2]-$matches-$nummatches+1;
			}
			else{
			    #Start is alignment start (s1) + matching columns 
			    $start = $aln->[1]+$matches+$nummatches-1;
			}
		    }
		    print "START genomic=$start\n" if($DEBUG > 1);
		    die if($start<$aln->[1]);
		    die if($start>$aln->[2]);
		    $startcount = $currcount+($columnstart-$currcount);
		    $foundstart=1;
		}
		if($count+$currcount>=$columnend){
		    if(!$foundend){
			print "FOUNDEND $currcount $matches $count$char\n" if($DEBUG > 1);
			if($char eq 'M'){
			    print "END IN MATCH $orient\n" if($DEBUG > 1);
			    if($orient eq '-'){
				$end = $aln->[2]-$matches-($columnend-$currcount);
			    }
			    else{
				$end = $aln->[1]+$matches+($columnend-$currcount);
			    }
			}
			else{
			    #Report last matching position
			    print "END INGAP $orient\n" if($DEBUG > 1);
			    #Report last matching position
			    #Last matching position is defined by last overlapping M interval between 
			    #query and current sequence
			    my $intersectbv = new Bit::Vector($querybv->Size());
			    $intersectbv->Intersection($querybv,$currbv);
			    die if($currbv->bit_test($currcount+1));
			    die if($intersectbv->bit_test($currcount+1));
			    my($imin,$imax) = $intersectbv->Interval_Scan_dec($currcount+1);
			    my($cmin,$cmax) = $currbv->Interval_Scan_dec($currcount+1);
			    my $nummatches=0;
			    if(! defined $cmax || !defined $imax){
				die "INGAP Can't find matching position in $aln->[0] < columnend:$columnend. No last matching position\n" if($DEBUG > 1);
			    }
			    else{
				for(my $i=$cmax;$i>=$imax;$i--){
				    if($currbv->bit_test($i)){
					$nummatches++;
				    }
				}
			    }
			    print "INGAP $nummatches until next match between query and $aln->[0] at column $cmax-$imax\n" if($DEBUG > 1);
			    if($orient eq '-'){
				$end = $aln->[2]-$matches;
			    }
			    else{
				$end = $aln->[1]+$matches;
			    }
			}
			print "END genomic=$end\n" if($DEBUG > 1);
			die if($end>$aln->[2]);
			die if($end<$aln->[1]);
			$foundend=1;
			last;
		    }
		}
	    }
	    if($char eq 'M'){
		$matches+=$count;
	    }
	    $currcount+=$count;
	}
	die "Could not find start or end " if(!$foundstart || !$foundend);
    }
    else{
	if($orient eq '-'){
	    $start=$aln->[2]-$columnend-1;
	    $end=$aln->[2]-$columnstart;
	}
	else{
	    $start=$aln->[1]+$columnstart-1;
	    $end=$aln->[1]+$columnend;
	}
    }
    #returning fmin<fmax
    if($orient eq '-'){
	($start,$end) = ($end,$start);
    }
     return ($start,$end,$currbv);
}

sub get_cigs {
    my($cig) = @_;
    my @chars = split /\d+/,$cig;
    my @counts = split /[MXIDG]/,$cig;

    my @cigs;
    my $columncount;
    for(my $i=0;$i<@counts;$i++){
	my $c = $chars[$i+1];
	die "Invalid cigar str $c\n" if($c !~ /[MXIDG]/);
	my $currcount = $counts[$i];
	die "Invalid count $currcount\n" if($currcount !~ /\d+/);
	push @cigs,[$currcount,$c];
	$columncount+=$currcount;
    }
    return (\@cigs,$columncount);
}

#return alnobj,bitvector,width
sub getAlignment{
    my($self,$name) = @_;
    die "Bad alignment $name" if(ref $name);
    if(exists $self->{_alignments}->{$name}){
	return @{$self->{_alignments}->{$name}};
    }
    else{
	print "Alignment $name not found\n";
	return undef;
    }
}

#Assumes one interval per genome, per alignment
sub getAlignedInterval{
    my($self,$align_name,$seqname) = @_;
    my $alignobj = $self->{_alignments}->{$align_name}->[0];
    foreach my $alni (@$alignobj){
	if($alni->[0] eq $seqname){
	    return $alni;
	}
    }
    print "#Can't find $seqname on alignment $align_name\n";
    return undef;
}


#
#Returns closed interval [$startcol,$endcol]
sub getAlignmentMatrix {
    my($self,$align_name,$startcol,$endcol,$db) = @_;
    if(!$startcol){
	$startcol=0;
    }
    
    die "Can't find alignment $align_name" if(!exists $self->{_alignments}->{$align_name});    
    my ($alignobj,$gv,$align_width) = @{$self->{_alignments}->{$align_name}};
    die "Bad input columns $startcol-$endcol $startcol >$align_width || $endcol > $align_width" if($startcol >$align_width || $endcol > $align_width);

    #populate alignment matrix
    my $row=0;
    my $matrix=[];
    my $seqmatrix=[];
    my @names;

    my $alni;

    if(!$endcol){
	$endcol=$align_width;
    }

    foreach my $alni (@$alignobj){
	my $matchcount=0;
	my $column=1;
	my $mstr = '-'x ($align_width+1);
	my $sstr = '-'x ($align_width+1);
	push @names,$alni->[0];
	my ($cigs,$columncount) = &get_cigs($alni->[4]);

	foreach my $c (@$cigs){
	    my($count,$char) = @$c;
	    if($char eq 'M'){
		my $mmstr = '.' x $count;
		die if(length($mmstr)!=$count);
		substr($mstr,$column) = $mmstr;
		my $seqobj = $db->get_Seq_by_id($alni->[0]);
		if($seqobj){
		    die if($alni->[1]>$alni->[2]);
		    my $str;
		    #print "$alni->[1]+$matchcount,$alni->[1]+$matchcount+$count-1\n";
		    if($alni->[3] eq '+'){
			$str = $seqobj->subseq($alni->[1]+$matchcount+1,$alni->[1]+$matchcount+$count);
		    }
		    else{
			#Note cigar always denotes offset from alignment start
			#In '-' orient, cigar starts from $alni->[2]--->$alni->[1]
			$str = revcom($seqobj->subseq($alni->[2]-$matchcount-$count+1,$alni->[2]-$matchcount))->seq();
		    }
		    die length($str)." != $count" if(length($str)!=$count);
		    die if(length($str)!=length($mmstr));
		    substr($sstr,$column,length($str)) = $str;
		    die if(substr($sstr,$column,length($str)) ne $str);
		}
		else{
		    die "Can't find seq $alni->[0]\n";
		}
		$matchcount+=$count;
	    }
	    else{
		my $mmstr = '-' x $count;
		die if(length($mmstr)!=$count);
		substr($mstr,$column) = $mmstr;
	    }
	    $column+=$count;
	}
	$matrix->[$row]=$mstr;
	$seqmatrix->[$row]=$sstr;
	$row++;
    }

    die "Invalid range $startcol-$endcol" if($endcol < $startcol);

    my $retmatrix=[];
    my $retseqmatrix=[];
    for(my $i=0;$i<@$matrix;++$i){ 
	$retmatrix->[$i] = substr($matrix->[$i],$startcol,$endcol-$startcol+1);
	$retseqmatrix->[$i] = substr($seqmatrix->[$i],$startcol,$endcol-$startcol+1);
	die "Bad sequence $retmatrix->[$i]" if(length ($retmatrix->[$i])<1);
	die "Bad sequence $retseqmatrix->[$i]" if(length ($retseqmatrix->[$i])<1);
    }
    #remove same characters
    #this is really going to really slow in perl this wa
    for(my $i=1;$i<@$retmatrix;++$i){ 
	for(my $j=0;$j<length($retseqmatrix->[$i]);$j++){
	    my $topchar = substr($retseqmatrix->[0],$j,1);
	    if(substr($retseqmatrix->[$i],$j,1) ne $topchar){
		substr($retmatrix->[$i],$j,1) = substr($retseqmatrix->[$i],$j,1);
	    }
	    else{
		if($topchar eq '-'){
		    substr($retmatrix->[$i],$j,1) = '-';
		}
		else{
		    substr($retmatrix->[$i],$j,1) = '.';
		}
	    }
	}
    }
    return ($retmatrix,$retseqmatrix,\@names);
}

sub contains{
    my($self,$align_name,$qseqname,$coord1,$coord2) = @_;
    die "Can't find alignment $align_name" if(!$align_name || !exists $self->{_alignments}->{$align_name});
    my $alignobj = $self->{_alignments}->{$align_name}->[0];
    foreach my $alni (@$alignobj){
	if($alni->[0] eq $qseqname){
	    if($coord1<$alni->[1] || $coord1>$alni->[2]){
		#print "Start position $coord1 is not contained in interval $qseqname:$alni->[1]-$alni->[2]\n";
		return 0;
	    }
	    elsif($coord2<$alni->[1] || $coord2>$alni->[2]){
		#print "End position $coord2 is not contained in interval $qseqname:$alni->[1]-$alni->[2]\n";
		return 0;
	    }
	    return 1;
	}
    }
    
}

#mappedfeats in the form [name,seq,start,end]
sub printAlignment{
    my($self,$aln,$startcol,$endcol,$db,$mappedfeats) = @_;
    die "Must specify Bioperl database $db that contains sequence data" if(!$db);
    die "Must specify startcol, endcol $startcol-$endcol" if(!$startcol || !$endcol);
    my($alignobj,$alignment_bv,$align_width) = @{$self->{_alignments}->{$aln}};

    #$mmatrix,$seqmatrix are relative to $startcol, index starting at 0
    my ($mmatrix,$seqmatrix,$names) = $self->getAlignmentMatrix($aln,$startcol,$endcol,$db);
    my $COL_WIDTH=100;
    my $atree = new AlignmentTree();
    my $features = {};
    foreach my $feat (@$mappedfeats){
	$atree->insert(@$feat);
	$features->{$feat->[1]} = [$feat->[0]->[0]->[1],$feat->[0]->[0]->[2]];
    }
    for(my $j=0;$j<=(($endcol-$startcol)/$COL_WIDTH);$j++){
	my $s=$j*$COL_WIDTH+1;
	my $e=$s+$COL_WIDTH-1;
	$e = ($e>($endcol-$startcol+1)) ? ($endcol-$startcol+1) : $e;
	my @coords;
	#offset into full alignment $aln
	my $absstartcol = $s+$startcol-1;
	my $absendcol = $e+$startcol-1;
	for(my $i=0;$i<@$names;$i++){
	    my($alni) = $self->getAlignedInterval($aln,$names->[$i]);
	    my($start,$end) = &columntocoords($alni,$absstartcol,$absendcol);
	    push @coords,[$start,$end,$alni->[3]];
	    ($start,$end) = ($alni->[3] eq '-') ? ($end,$start) : ($start,$end);
	    my $displaystr;
	    if($i==0){
		$displaystr = substr($seqmatrix->[$i],$s-1,$e-$s+1);
	    }
	    else{
		$displaystr = substr($mmatrix->[$i],$s-1,$e-$s+1);
	    }
	    if($self->{debug}){
		printf("%30.30s %7s %11s %-30s %7s %11s\n",
		       "$names->[$i]:$alni->[3]",
		       $start,
		       "col:$absstartcol",
		       $displaystr,
		       $end,
		       "col:$absendcol");
	    }
	    else{
		printf("%30.30s %7s %11s %-30s %7s %11s\n",
		       "$names->[$i]:$alni->[3]",
		       $start,
		       "",
		       $displaystr,
		       $end,
		       "");
	    }
	}
	printf("%".$COL_WIDTH.".".$COL_WIDTH."s","$aln col:$absstartcol-$absendcol\n");
	for(my $i=0;$i<@$names;$i++){
	    #Show all matching features that intersect $coords[$i]->[0],$coords[$i]->[1]
	    if($coords[$i]->[1]-$coords[$i]->[0]>1){
		my @res = $atree->intersect($names->[$i],$coords[$i]->[0],$coords[$i]->[1],'gene');
		foreach my $r (@res){
		    
		    #offset into full alignment $aln
		    my($cs,$ce) = &coordstocolumn($alignobj,$names->[$i],$r->[2],$r->[3]);
		    #print "$r->[0] coords:$r->[2],$r->[3] $cs,$ce $cs-$absstartcol $absendcol-$ce\n";
		    my $leadinggap = 'X'x($cs-$absstartcol);
		    my $trailinggap = 'X'x($absendcol-$ce);
		    my $displaystr;
		    my $startcodonstr;
		    my $stopcodonstr;
		    my $displaytoken;
		    
		    #TODO, REFACTOR into a matrix. this impl doesn't support condons that span row bounds
		    #currently only viz start,stop codons at beginning/end of alignment
		    #print "$r->[2] <= $features->{$r->[0]}->[0] && $r->[3] >= $features->{$r->[0]}->[0]\n";
		    if($r->[2] <= $features->{$r->[0]}->[0] && $r->[3] >= $features->{$r->[0]}->[0]){
			if($coords[$i]->[2] eq '-'){
			    if($r->[6] eq '-'){
				$stopcodonstr = 'TAA';#substr($seqmatrix->[$i],$cs-$absstartcol-($COL_WIDTH-3)+($j*$COL_WIDTH),3);
				$displaytoken .= 'STOP1<--';
			    }
			    else{
				$stopcodonstr = 'CAT';#substr($seqmatrix->[$i],$cs-$absstartcol-($COL_WIDTH-3),3);
				$displaytoken .= '<--START1';
			    }
			}
			else{
			    if($r->[6] eq '-'){
				$startcodonstr = 'TTA';#substr($seqmatrix->[$i],$cs-$absstartcol+($j*$COL_WIDTH),3);
				$displaytoken  .= 'STOP2<--';
			    }
			    else{
				$startcodonstr = 'ATG';#substr($seqmatrix->[$i],$cs-$absstartcol+($j*$COL_WIDTH),3);
				$displaytoken .= 'START2-->';
			}
			}
		    #TODO trim to row
		    }
		    if($r->[2] <= $features->{$r->[0]}->[1] && $r->[3] >= $features->{$r->[0]}->[1]){
			if($coords[$i]->[2] eq '-'){
			    #$stopcodonstr = substr($seqmatrix->[$i],$absendcol-$ce-($COL_WIDTH-3),3);
			    if($r->[6] eq '-'){
				#$startcodonstr = substr($seqmatrix->[$i],$absendcol-$ce-($COL_WIDTH-3),3);
				$startcodonstr = 'ATG';#substr($seqmatrix->[$i],$cs-$absstartcol+($j*$COL_WIDTH),3);
				$displaytoken .= 'START3-->';
			    }
			    else{
				$startcodonstr = 'TTA';#substr($seqmatrix->[$i],$absendcol-$ce-($COL_WIDTH-3),3);
				$displaytoken .= 'STOP3<--';
			    }
			}
			else{
			    if($r->[6] eq '-'){
				$stopcodonstr = 'CAT';#substr($seqmatrix->[$i],$absendcol-$ce-3,3);
				$displaytoken .= '<--START4';
			    }
			    else{
				$stopcodonstr = 'TAA';#substr($seqmatrix->[$i],$absendcol-$ce-3,3);
				$displaytoken .= 'STOP4<--';
			    }
			}
			#TODO trim to row
		    }
		    my $spacer = '_'x($ce-$cs+1-length($startcodonstr)-length($stopcodonstr));
		    
		    $displaystr =  $startcodonstr.$spacer.$stopcodonstr;#;substr($seqmatrix->[$i],$cs-1,$ce-$cs+1);
		    die if(length($displaystr) > $COL_WIDTH);
		    my ($feat_start,$feat_end) = ($r->[2],$r->[3]);
		    
		    #TODO determine frame and print ~ only every 3 codons in frame if in frame
		    #otherwise
		    my $frame;
		    if($coords[$i]->[2] eq '-'){
			if($r->[6] eq '-'){
			    die if($features->{$r->[0]}->[1] < $feat_end);
			    $frame = (($features->{$r->[0]}->[1] - $feat_end)%3);
			}
			else{
			    die if($features->{$r->[0]}->[1] < $feat_end);
			    $frame = (($features->{$r->[0]}->[1] - $feat_end)%3);
			}
		    }
		    else{
			if($r->[6] eq '-'){
			    die "$feat_start < $features->{$r->[0]}->[0] $r->[0]" if($feat_start < $features->{$r->[0]}->[0]);
			    $frame = (($feat_start - $features->{$r->[0]}->[0])%3);
			}
			else{
			    die if($feat_start < $features->{$r->[0]}->[0]);
			    $frame = (($feat_start - $features->{$r->[0]}->[0])%3);
			}
		    }
		    ($feat_start,$feat_end) = ($coords[$i]->[2] eq '-') ? ($feat_end,$feat_start) : ($feat_start,$feat_end);


		    my $m=$frame;
		    my $fulldisplaystr = $leadinggap.$displaystr.$trailinggap;
		    for(my $k=length($leadinggap);$k<length($leadinggap)+length($displaystr);$k++){
			my $idx=$k;
			if(substr($mmatrix->[$i],$s-1+$k,1) ne '-'){
			    $m++;
			    if($m%3==0){
				#Don't overwrite start,stop codons
				#if($idx>length($startcodonstr)+length($leadinggap)
				#   & $idx<(length($leadinggap)+length($displaystr)-length($stopcodonstr))){
				substr($fulldisplaystr,$idx,1) = '|';#$frame;
				#}
			    }
			    else{
				#substr($fulldisplaystr,$idx,1) = substr($mmatrix->[$i],$s-1+$k,1);	
			    }
			}
			else{
			    if(substr($mmatrix->[$i],$s-1+$k,1) eq '-'){
				substr($fulldisplaystr,$idx,1) = substr($mmatrix->[$i],$s-1+$k,1);	
			    }
			}
		    }		    
		    die if($r->[6] ne $r->[7]);
		    printf("%30.30s %7s %11s %-30s %7s %11s\n",
			   $r->[0].":$r->[6]",
			   $feat_start,
			   $displaytoken,
			   $fulldisplaystr,
			   $feat_end,
			   $displaytoken);
		}
	    }
	}
	printf("%".$COL_WIDTH.".".$COL_WIDTH."s","ANNOTATIONS\n");
    }
}

sub printAlignmentDebug{
    my($alignobj,$handle) = @_;
    foreach my $alni2 (@$alignobj){
	if(!$handle){
	    $handle=\*STDOUT;
	}
	print $handle "#ALIGNOBJ $alignobj ",join(' ',@$alni2),"\n";
    }
}

sub removeOverlaps{
    my($self,$alignments,$qseqname) = @_;
    my @alns;
    my @results;
    my %contained;
    my %overlaps;
    foreach my $align_name (@$alignments){
	my $alni = $self->getAlignedInterval($align_name,$qseqname);
	if($align_name =~ /$aligntoken/){
	    push @alns,[$align_name,$alni->[1],$alni->[2],$alni->[2]-$alni->[1]];
	}
    }
    my @sortedalns = sort {$b->[3] <=> $a->[3]} @alns;
    for(my $i=0;$i<@sortedalns;$i++){
	my $ifmin = $sortedalns[$i]->[1];
	my $ifmax = $sortedalns[$i]->[2];
	for(my $j=$i+1;$j<@sortedalns;$j++){
	    my $jfmin = $sortedalns[$j]->[1];
	    my $jfmax = $sortedalns[$j]->[2];
	    if($jfmin>=$ifmin && $jfmax <=$ifmax){
		print "Marking $sortedalns[$j]->[0] contained $jfmin>=$ifmin && $jfmax <=$ifmax in $sortedalns[$i]->[0]\n" if($DEBUG);
		$contained{$j}++;
	    }
	    else{
		if($jfmin>=$ifmin && $jfmin <=$ifmax){
		    $overlaps{$j}++;
		}
		if($jfmax>=$ifmin && $jfmax <=$ifmax){
		    $overlaps{$j}++;
		}
	    }
	}
    }
    if(scalar(keys %overlaps)>0){
	print "#WARNING removing some alignments with overlaps\n" if($DEBUG);;
    }
    if(scalar(keys %contained)>0){
	print "#WARNING removing some alignments that are fully contained\n" if($DEBUG);;
	for(my $i=0;$i<@sortedalns;$i++){
	    if(!exists $contained{$i}){
		push @results,$sortedalns[$i]->[0];
	    }
	    else{
		print "#WARNING removing contained alignment $sortedalns[$i]->[0]\n" if($DEBUG);;
	    }
	}
	return @results;
    }
    else{
	return @$alignments;
    }
}
1;
