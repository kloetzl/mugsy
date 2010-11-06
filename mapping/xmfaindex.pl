#!/usr/bin/perl
#
#./mafindex.pl mugsyindex < mugsy.out
#Adds an MAF formatted file to a MUGSY formatted index
#Each alignment is saved as type 'alignment'
#
use strict;
use lib '/usr/local/projects/angiuoli/developer/sangiuoli/mugsy/trunk';
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

my $index=0;
my $seqlookup = {};
if(-e $ARGV[1]){
    open FILE,"$ARGV[1]" or die "Can't open file $ARGV[1]";
    while(my $line=<FILE>){
	my($seq) = ($line =~ /\>?(\S+)/);
	$seqlookup->{++$index} = $seq;
    }
    close FILE;
}

my $currscore;
my $block = [];
my $k=0;
my $label=0;

my $seqname;
my $start;
my $end;
my $orient;
my @seqinfo;
while(my $line=<STDIN>){
    if($line =~ /^=/){
	if(defined $seqname && $start>0){
	    my ($cigar,$len) = &get_cigar(join('',@seqinfo));
	    die "Bad match length $len in cigar $cigar" if ($end-$start+1 != $len);
	    #Convert alignment to zero start, interbase coordinates
	    push @$block,[$seqname,$start-1,$end,$orient,$cigar];
	    print "Adding aligned sequence $seqname $start-1,$end,$orient to alignment MAUVE_$label\n";
	}
	$atree->insert($block,"MAUVE_$label","alignment") if(scalar(@$block));
	$label++;
	$block=[];
	$seqname=undef;
	@seqinfo=();
    }
    elsif($line =~ /^>\s+(\d+)\:(\d+)-(\d+)\s+([\+\-])\s+(\S+)/){
	if(defined $seqname && $start>0){
	    my ($cigar,$len) = &get_cigar(join('',@seqinfo));
	    die "Bad match length $len in cigar $cigar" if ($end-$start+1 != $len);
	    push @$block,[$seqname,$start-1,$end,$orient,$cigar];
	    print "Adding aligned sequence $seqname $start-1,$end,$orient to alignment MAUVE_$label\n";
	}
	my $seqid = $1;
	$start = $2;
	if($start>0){
	    $end = $3;
	    #XMFA format start always < end 
	    die "Invalid coordinates $start-$end" if($start>$end);
	    #Relative orientation of the alignment
	    $orient = $4;
	    my $file = $5;
	    $seqname = $file;
	    if(exists $seqlookup->{$seqid}){
		$seqname = $seqlookup->{$seqid};
	    }
	    else{
		#Hack for strep pneumo xmfa files
		$seqname =~ s/\.fsa//g;
	    }
	}
	@seqinfo=();
    }
    else{
	if(defined $seqname){
	    chomp $line;
	    push @seqinfo,$line;
	}
    }
}
$atree->insert($block,"MAUVE_$label","alignment") if(scalar(@$block));
print STDERR "Writing index to $ARGV[0]\n";
$atree->serialize($ARGV[0]);

sub get_cigar{
    my($seqs) = @_;
    my $cig;
    my $len=0;
    my $mlen=0;
    my @chars = split(//,$seqs);
    my $count=0;
    my $curr=0; #1 - match, 2 - gap
    foreach my $c (@chars){
	#match char
	if($c ne '-'){
	    if($curr==2){
		#in gap
		#write prev gap
		$cig .= $count."X";
		$count=0;
	    }
	    #in match
	    $count++;
	    $curr=1;
	}
	else{
	    #gap char
	    if($curr==1){
		#in match
		#write prev gap
		$cig .= $count."M";
		$len += $count;
		$count=0;
	    }
	    #in gap
	    $count++;
	    $curr=2;
	}
    } 
    if($curr==1){
	#in gap
	#write prev gap
	$cig .= $count."M";
	$len += $count;
    }
    if($curr==2){
	#in gap
	#write prev gap
	$cig .= $count."X";
    }
    return ($cig,$len);
}
