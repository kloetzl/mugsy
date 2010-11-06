#!/usr/bin/perl
#./reportvariants.pl index fasta

use strict;
use Bio::Perl;
use Bio::DB::Fasta;
use Bio::Seq;
use lib '/usr/local/projects/angiuoli/developer/sangiuoli/mugsy/trunk/mapping/';
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use AlignmentTree;

my %options;
my $results = GetOptions (\%options, 
			  'gap_window|g=s',
			  'display_window|d=s',
			  'gaps_allowed|a=s') || pod2usage(-verbose => 1);

pod2usage(-verbose=>1) if($options{'help'});

my $atree = AlignmentTree::deserialize($ARGV[0]);

my $db = Bio::DB::Fasta->new($ARGV[1],'-reindex'=>1); 

my $gapthreshold=0;
if(exists $options{'gaps_allowed'}){
    $gapthreshold = $options{'gaps_allowed'};
}
my $gap_window=5;
if(exists $options{'gap_window'}){
    $gap_window = $options{'gap_window'};
}
my $display_window=5;
if(exists $options{'display_window'}){
    $display_window = $options{'display_window'};
}


shift @ARGV;
shift @ARGV;

my $pwseqs = {};
my $refname = shift @ARGV;
foreach my $seq (@ARGV){
    $pwseqs->{$seq}++;
}

open VFILE,"+>$$.pwvariants.out" or die "Can't open file pwvariants.out";
open SFILE,"+>$$.snpvariants.out" or die "Can't open file snpvariants.out";
foreach my $alnname (sort {$a cmp $b} keys %{$atree->{_alignments}}){
    my($alnobj,$aln_bv,$align_width) = @{$atree->{_alignments}->{$alnname}};
    my ($mmatrix,$seqmatrix,$names) = $atree->getAlignmentMatrix($alnname,1,$align_width,$db);
    if(@$seqmatrix > 1){
	#print STDERR "Checking alignment $alnname $align_width ",scalar(@$seqmatrix),"\n";
	
	my $ngaps;
	my $nmismatches;
	my $variants = {};
	my $seqvariants = {};
	my $refidx;
	for(my $i=0;$i<@$seqmatrix;$i++){
	    if($names->[$i] eq $refname){	
		$refidx=$i;
	    }
	}
#Matrix cols start at 0
	for(my $j=0;$j<$align_width;$j++){
	    my $b;
	    my $refbp = lc(substr($seqmatrix->[$refidx],$j,1));
	    for(my $i=0;$i<@$seqmatrix;$i++){
		if($i ne $refidx){
		    my $currbp = lc(substr($seqmatrix->[$i],$j,1));
		    if($currbp ne $refbp && $currbp !~ /[yskrmwnw]/){
			$variants->{$j}++;
			$seqvariants->{$i}->{$j}++;
		    }
		}
		#print "$b=$currbp " if($b ne '-' && $currbp ne '-');
	    }
	}
	#print STDERR "variants ",scalar(keys %$variants),"\n";
	foreach my $col (sort {$a <=> $b} keys %$variants){
	    my $gaps=0;
	    for(my $i=0;$i<@$seqmatrix;$i++){
		my $start = $col - $gap_window;
		$start = 0 if($start < 0);
		my $end = $col + $gap_window;
		$end = $align_width if($end > $align_width);
		$gaps+= (substr($seqmatrix->[$i],$start,$end-$start+1) =~ tr/\-/\-/);
	    }
	    if($gaps<=$gapthreshold){
		my $refc;
		for(my $i=0;$i<@$seqmatrix;$i++){
		    my $start = $col - $display_window;
		    $start = 0 if($start < 0);
		    my $end = $col + $display_window;
		    $end = $align_width if($end > $align_width);
		    my($alni) = $atree->getAlignedInterval($alnname,$names->[$i]);
		    my $colstart = 1+$start;
		    my $colend = $colstart;
		    my($startc,$endc) = AlignmentTree::columntocoords($alni,$col+1,$col+1);
		    $refc = $startc if($names->[$i] eq "$refname");
		    #AlignmentTree::printAlignmentDebug($alnobj);
		    printf("%10s %s\tcoords:%d-%d\n",$names->[$i],lc(substr($seqmatrix->[$i],$start,$end-$start+1)),$startc,$endc);
#, substr($seqmatrix->[$i],$start,$end-$start),"\n";
		    
		    if($names->[0] eq "$refname" && exists $pwseqs->{$names->[$i]} && $seqvariants->{$i}->{$col}){
			print SFILE "$names->[$i]\t$refname\t$refc\t",$refc+1,"\t",uc(substr($seqmatrix->[0],$col,1)),"\n";
			print VFILE "$names->[$i]\t$refc\t",$refc+1,"\t",substr($seqmatrix->[0],$col,1),"/",substr($seqmatrix->[$i],$col,1),"\t$names->[$i]\t$startc-$endc\n";
		    }
		}
		printf("%10s      ^     \n");
		print "\n";
	    }
	}
    }
}
close VFILE;
close SFILE;
