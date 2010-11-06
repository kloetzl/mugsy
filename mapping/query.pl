use strict;
use AlignmentTree;
use Storable qw(store retrieve);
use Data::Dumper;

$Storable::Deparse = 1;
$Storable::Eval = 1;

my $atree;
if(-e $ARGV[0]){
    $atree = retrieve($ARGV[0]);
}

my @results = $atree->intersect($ARGV[1],$ARGV[2],$ARGV[3]);

foreach my $r (@results){
    print "INTERSECT RESULT ",join(' ',@$r),"\n";
}

