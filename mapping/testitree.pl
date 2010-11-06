#!/usr/bin/perl

use strict;
use IntervalTree;
use AlignmentTree;
use Data::Dumper qw(Dumper);

#remove only using for revcom
use Bio::Perl;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::Tools::CodonTable;

#Assumptions fmin<fmax,colstart<colend
#Genome coordinate system is 0 start, interbase coordinates. Feature length = fmax-fmin
#Alignment coordinate system is 1 start counting bases. Feature length is fmax-fmin+1.

#Test cases
#Sequence 1 ...AATTGGCCAA...
#Sequence 2 ...AATTGGCCAA...
#Sequence 3 ...AATTGGCCAA...

#Alignment 1 S1,S2,S3 +,+,+
#Alignment 2 +,-,+
#Alignment 3 -,-,+

#Test feature1 orient='+' fmin=102 fmax=107 'TTGGC'
#Test feature2 orient='-' fmin=103 fmax=108 'GCCAA'
#
#+ Alignment, + annotation end5<end3 colorient '+' fmin -> coords increasing -> fmax
#Eg. feature1
#100 1 AATTGGCCAA 10 110 
#100 1 AATTGGCCAA 10 110 
#        TTGGC    
#         GCCAA
#col   123456789
#query:fmin=102,fmax=107 strand +
#result:colstart=3,coldend=7,revcomp=0

#- Alignment, - annotation end3<end5 colorient '+' fmax -> coords decreasing -> fmin
#Eg. feature2
#110 1 TTGGCCAATT 10 100 
#110 1 TTGGCCAATT 10 100 
#         GCCAA
#         CGGTT         
#col   123456789
#query:fmin=102,fmax=107 strand -
#result:colstart=4,colend=8,revcomp=0

#+ Alignment, - annotation end3<end5 colorient '-' fmin -> coords increasing -> fmax. revcom matching interval
#Eg. feature2
#100 1 AATTGGCCAA 10 110 
#100 1 AATTGGCCAA 10 110 
#        AACCG - reversed 107-102
#col   123456789
#query:fmin=102,fmax=107 strand - 
#result:colstart=3,colend=7,revcomp=1

#- Alignment, + annotation end5<end3 colorient '-' fmax -> coords decreasing -> fmin. revcom matching interval
#120 1 TTGGCCAATT 10 100 
#110 1 TTGGCCAATT 10 100 
#         CGGTT - reversed 107-102
#col   123456789
#query:fmin=102,fmax=107 strand + 
#result:colstart=4,colend=8,revcomp=1



my @alignments = ([
		   ['genome1',10,1000,'+','900M100X','g1'],
		   ['genome2',100,900,'+','100X800M100X','g2'],
		   ['genome3',350,1350,'+','1000M','g3'],
		   ],
		  [
		   ['genome1',20,2000,'+','1820M180X','g1'],
		   ['genome2',200,900,'+','180X700M1120X','g2'],
		   ['genome3',450,2350,'+','100X1900M','g3'],
		   ['genome4',450,2350,'+','100X1900M','g4']
		   ]
		  );

my @alignqueries = (["genome1",1010,1020],
		    ["genome2",500,720]
		    );

my @intervals = ([10,1000,1,'+'],
		 [100,900,2,'+'],
		 [350,10000,3,'+']);

my @intqueries = ([1010,1020],
	       [500,720]
	       );

my @filter = ('g1','g4');

#
#Test intervaltree
my $tree = new IntervalTree(1,1000000);
foreach my $i (@intervals){
    $tree->insert(@$i);
}

#for(my $i=10000;$i>=0;$i--){
#    print "$i\n" if($i%1000==0);
#    $tree->insert($i,$i+1,$i);
    #print Dumper($tree),"\n";
#}

#for(my $i=0;$i<10000;$i++){
#    print "$i\n" if($i%1000==0);
#    $tree->insert($i,$i+1,$i);
#}

foreach my $q (@intqueries){
    print "QUERY ",join(' ',@$q),"\n";
    my @results = $tree->intersect(@$q);
    foreach my $r (@results){
	print "RESULT $r\n";
    }
}

#
#Test alignment tree
my $atree = new AlignmentTree();

my $k=0;
foreach my $a (@alignments){
    $atree->insert($a,"MAUVE$k","MAUVE");
    $k++;
}
print "Alignmenttree intersect queries\n";
foreach my $q (@alignqueries){
    print "QUERY ",join(' ',@$q),"\n";
    my @results = $atree->intersect(@$q);
    
    foreach my $r (@results){
	print "INTERSECT RESULT ",join(' ',@$r),"\n";
    }
    print "DONE\n";
}
print "Alignmenttree map()\n";
foreach my $q (@alignqueries){
    print "QUERY ",join(' ',@$q),"\n";
    my @results = $atree->map(@$q);
    
    foreach my $r (@results){
	print "MAP RESULT ",join(' ',@$r),"\n";
    }
}
print "DONE\n";

print "Adding filter ",join(',',@filter),"\n";
$atree->filter(@filter);

foreach my $q (@alignqueries){
    print "QUERY ",join(' ',@$q),"\n";
    my @results = $atree->intersect(@$q);
    
    foreach my $r (@results){
	print "INTERSECT RESULT ",join(' ',@$r),"\n";
    }
}


#TEST 1 +,+ alignment mapped features on opposing strands
#
#Test alignment tree

open FILE,">/tmp/$$.testing" or die "Can't open file /tmp/$$.testing";
print FILE <<_FASTAEND;
>genome1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
AATTGGCCAANNNN
>genome2
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
AATTGGCCAANNNN
_FASTAEND

    ;
close FILE;

my $db = Bio::DB::Fasta->new("/tmp/$$.testing",'-reindex'=>1); 

my $atree = new AlignmentTree();

#0 1 2 3 4 5 6 7 8 9 10
# A A T T G G C C A A

my @alignments2 = ([
		    ['genome1',100,110,'+','10M','g1'], #AATTGGCCAA
		    ['genome2',100,110,'+','10M','g2']  #AAATTGGCCA
		    ],
		   [
		    ['genome1',100,110,'+','5M1X5M','g1'], #AATTGGCCAA
		    ['genome2',100,110,'+','5M1X5M','g2']  #AAATTGGCCA
		    ]);

my $expectedalns = ['AATTGGCCAA','AATTGGCCAA','AATTGGCCAA','AATTGGCCAA'];

my @features = ([['genome1',102,107,'+','5M']], #ATTGG
		[['genome2',102,107,'-','5M']]  #CCAAT
		);
my $expectedfeats = ['TTGGC','GCCAA','TTGGC','GCCAA'];


my $k=0;
my $alnidx=0;
foreach my $a (@alignments2){
    foreach my $f (@$a){
	#convert from 0 base to 1 base
	my $queryseq = $db->get_Seq_by_id($f->[0]);
	my $queryseqsubstr = $queryseq->subseq($f->[1]+1,$f->[2]);
	if($f->[3] eq '-'){
	    $queryseqsubstr = revcom($queryseqsubstr)->seq();
	}
	else{
	    $queryseqsubstr = $queryseq->subseq($f->[1]+1,$f->[2]);
	}
	if($queryseqsubstr =~ /N/){
	    die "ERROR unexpected alignment sequence $k found $queryseqsubstr =~ /N/\n";
	}
	if($queryseqsubstr ne $expectedalns->[$k]){
	    die "ERROR unexpected alignment sequence $k found $queryseqsubstr ne $expectedalns->[$k]\n";
	}
	$k++;
    }
    $atree->insert($a,"MAUVE$alnidx","MAUVE");
    $alnidx++;
}
$k=0;
$alnidx=0;
foreach my $f (@features){
    $atree->insert($f,'gene:'.$k,'gene');
    #convert from 0 base to 1 base
    my $queryseq = $db->get_Seq_by_id($f->[0]->[0]);
    my $queryseqsubstr = $queryseq->subseq($f->[0]->[1]+1,$f->[0]->[2]);
    if($f->[0]->[3] eq '-'){
	print "REVCOM $queryseqsubstr\n";
	$queryseqsubstr = revcom($queryseqsubstr)->seq();
    }
    if($queryseqsubstr =~ /N/){
	die "ERROR unexpected sequence found $queryseqsubstr =~ /N/\n";
    }
    if($queryseqsubstr ne $expectedfeats->[$k]){
	die "ERROR unexpected sequence found $queryseqsubstr ne $expectedfeats->[$k]\n";
    }
    print "QUERYSEQ $f->[0]->[1]-$f->[0]->[2] $f->[0]->[3] ",$queryseqsubstr,"\n";
    $k++;
}

print "INTERSECT TEST1\n";
$k=0;
foreach my $f (@features){
    my @results = $atree->intersect($f->[0]->[0],$f->[0]->[1],$f->[0]->[2],'gene');
    die "More results than expected" if(scalar(@results)>1);
    my $r = $results[0];
    print "genome1 INTERSECT RESULT ",join(' ',@$r),"\n";

    my $queryseq = $db->get_Seq_by_id($r->[1]);
    my $queryseqsubstr = $queryseq->subseq($r->[2]+1,$r->[3]);
    if($r->[6] eq '-'){
	print "REVCOM $queryseqsubstr\n";
	$queryseqsubstr = revcom($queryseqsubstr)->seq();
    }
    if($queryseqsubstr =~ /N/){
	die "ERROR unexpected sequence found $queryseqsubstr =~ /N/ $r->[2]+1,$r->[3]\n";
    }
    if($queryseqsubstr ne $expectedfeats->[$k]){
	die "ERROR unexpected sequence found $queryseqsubstr ne $expectedfeats->[$k] $r->[2]+1,$r->[3]\n";
    }
    $k++;

    my @results1 = $atree->map($f->[0]->[0],$f->[0]->[1],$f->[0]->[2],'MAUVE');
    foreach my $r (@results1){
	print "MAP RESULT ",join(' ',@$r),"\n";
    }
}
$k=0;
my $alnidx=0;
#For all alignments
foreach my $a (@alignments2){
    my($alnobj,$bv,$width) = $atree->getAlignment("MAUVE$alnidx");
    my ($mmatrix,$seqmatrix,$names) = $atree->getAlignmentMatrix("MAUVE$alnidx",1,$width,$db);
    #See if we can map features into alignment matrix
    for(my $i=0;$i<@features;$i++){
	my $f = $features[$i];
	print "MATRIX MAUVE$alnidx ",join(',',@{$f->[0]}),"\n";
	my $flen = $f->[0]->[2]-$f->[0]->[1];
	die "Bad sequence $seqmatrix->[$i]" if(length ($seqmatrix->[$i])<1);
	#Returned seq length == input length + 1
	my($cs,$ce) = AlignmentTree::coordstocolumn($alnobj,$f->[0]->[0],$f->[0]->[1],$f->[0]->[2]);
	my $queryseqsubstr = substr($seqmatrix->[$i],$cs-1,$ce-$cs+1);
	$queryseqsubstr =~ s/\-//g;
	if($features[$i]->[0]->[3] eq '-'){
	    print "REVCOM $queryseqsubstr\n";
	    $queryseqsubstr = revcom($queryseqsubstr)->seq();
	}

	if($queryseqsubstr =~ /N/){
	    die "ERROR unexpected sequence found $queryseqsubstr =~ /N/\n";
	}
	if($queryseqsubstr ne $expectedfeats->[$k]){
	    die "ERROR unexpected sequence $k found $queryseqsubstr ne $expectedfeats->[$k]\n";
	}
	else{
	    print "Sequence $k $queryseqsubstr eq $expectedfeats->[$k] OK\n";
	}
	$k++;
    }
    $atree->printAlignment("MAUVE$alnidx",1,$width,$db);
	
    $alnidx++;
}
$alnidx=0;
foreach my $a (@alignments2){
    my($alnobj,$bv,$width) = $atree->getAlignment("MAUVE$alnidx");
    foreach my $f (@features){
	my @results1 = $atree->map($f->[0]->[0],$f->[0]->[1],$f->[0]->[2],'MAUVE');
	foreach my $r (@results1){
	    print "MAP RESULT ",join(' ',@$r),"\n";
	}
	$atree->printAlignment("MAUVE$alnidx",1,$width,$db,\@results1);
    }
    $alnidx++;
}
