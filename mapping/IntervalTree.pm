package IntervalTree;

#Adapted from bx-python

use strict;
use Math::Random qw(random_uniform);
use POSIX qw(ceil floor);
use Data::Dumper;

sub new{
    my $classname = shift;
    my $self = {};
    bless($self,$classname);
    my($start,$end,$name,$orient) = @_;
    $self->{'priority'} = ceil( (-1.0 / log(.5)) * log( -1.0 / (random_uniform(1,0,1) - 1)));
    $self->{'start'} = $start;
    $self->{'end'} = $end;
    $self->{'orient'} = $orient;
    $self->{'maxend'} = $self->{'end'};
    $self->{'minend'} = $self->{'end'};
    $self->{'left'} = undef;
    $self->{'right'} = undef;
    $self->{'name'} = $name;
    $self->{'default_func'} = sub {
	my $interval = shift; 
	return $interval->{'name'};
    };
    return $self;
}

sub insert{
    my($self,$start,$end,$name,$orient) = @_;
    die "Bad start-end $start-$end" if($end<$start);
    die "Bad orient $orient" if($orient ne '-' && $orient ne '+');
    my $root = $self;
    if($start > $self->{'start'}){
	if(defined $self->{'right'}){
	    $self->{'right'} = $self->{'right'}->insert($start,$end,$name,$orient);
	}
	else{
	    $self->{'right'} = new IntervalTree($start,$end,$name,$orient);
	}
	# rebalance tree
	if($self->{'priority'} < $self->{'right'}->{'priority'}){
	#    $root = $self->rotateleft();
	}
    }
    else{
	if(defined $self->{'left'}){
	    $self->{'left'} = $self->{'left'}->insert($start,$end,$name,$orient);
	}
	else{
	    $self->{'left'} = new IntervalTree($start, $end, $name,$orient);
	}
	# rebalance tree
	if($self->{'priority'} < $self->{'left'}->{'priority'}){ 
	#    $root = $self->rotateright();
	}
    }
    if(defined $root->{'right'} && defined $root->{'left'}){
	$root->{'maxend'} = ($root->{'end'}>$root->{'right'}->{'maxend'}) ? $root->{'end'} : $root->{'right'}->{'maxend'};
	$root->{'maxend'} = ($root->{'maxend'}>$root->{'left'}->{'maxend'}) ? $root->{'maxend'} : $root->{'left'}->{'maxend'};

	$root->{'minend'} = ($root->{'end'}<$root->{'right'}->{'minend'}) ? $root->{'end'} : $root->{'right'}->{'minend'};
	$root->{'minend'} = ($root->{'minend'}<$root->{'left'}->{'minend'}) ? $root->{'minend'} : $root->{'left'}->{'minend'};
    }
    elsif(defined $root->{'right'}){
	$root->{'maxend'} = ($root->{'end'}>$root->{'right'}->{'maxend'}) ? $root->{'end'} : $root->{'right'}->{'maxend'};
	$root->{'minend'} = ($root->{'end'}<$root->{'right'}->{'minend'}) ? $root->{'end'} : $root->{'right'}->{'minend'};
    }
    elsif(defined $root->{'left'}){
	$root->{'maxend'} = ($root->{'end'}>$root->{'left'}->{'maxend'}) ? $root->{'end'} : $root->{'left'}->{'maxend'};
	$root->{'minend'} = ($root->{'end'}<$root->{'left'}->{'minend'}) ? $root->{'end'} : $root->{'left'}->{'minend'};
    }
    return $root;
}

sub intersect{
    my($self,$start,$end,$func) = @_;
    die "$self->{'name'}:$end<$start not valid query" if($end<$start);
    my @results;

    $func = $self->{'default_func'} if(!$func);
    #print "CHECKING $start $end\n";
    if($start < $self->{'end'} && $end > $self->{'start'}){
	#print "Found $self->{'name'} $start <= $self->{'end'} && $end >= $self->{'start'}\n";
	push @results,$func->( $self );
    }
    if(defined $self->{'left'} && $start <= $self->{'left'}->{'maxend'}){
	push @results, $self->{'left'}->intersect( $start, $end, $func );
    }
    if(defined $self->{'right'} && $end >= $self->{'start'}){
	push @results, $self->{'right'}->intersect( $start, $end, $func );
    }

     return @results;
}

sub rotateright{
    my($self) = @_;
    die if(!exists $self->{'left'});
    die if(!exists $self->{'left'}->{'right'});
    my $root = $self;
    if(defined $self->{'left'}->{'right'}){
	$root = $self->{'left'};
	$self->{'left'} = $self->{'left'}->{'right'};
	$root->{'right'} = $self;
	if(defined $self->{'right'} && defined $self->{'left'}){
	    $self->{'maxend'} = ($self->{'end'}>$self->{'right'}->{'maxend'}) ? $self->{'end'} : $self->{'right'}->{'maxend'};
	    $self->{'maxend'} = ($self->{'maxend'}>$self->{'left'}->{'maxend'}) ? $self->{'maxend'} : $self->{'left'}->{'maxend'};
	    
	    $self->{'minend'} = ($self->{'end'}<$self->{'right'}->{'minend'}) ? $self->{'end'} : $self->{'right'}->{'minend'};
	    $self->{'minend'} = ($self->{'minend'}<$self->{'left'}->{'minend'}) ? $self->{'minend'} : $self->{'left'}->{'minend'};
	}
	elsif(defined $self->{'right'}){
	    $self->{'maxend'} = ($self->{'end'}>$self->{'right'}->{'maxend'}) ? $self->{'end'} : $self->{'right'}->{'maxend'};
	    $self->{'minend'} = ($self->{'end'}<$self->{'right'}->{'minend'}) ? $self->{'end'} : $self->{'right'}->{'minend'};
	}
	elsif(defined $self->{'left'}){ 
	    $self->{'maxend'} = ($self->{'end'}>$self->{'left'}->{'maxend'}) ? $self->{'end'} : $self->{'left'}->{'maxend'};
	    $self->{'minend'} = ($self->{'end'}<$self->{'left'}->{'minend'}) ? $self->{'end'} : $self->{'left'}->{'minend'};
	}
    }
    return $root;
}
sub rotateleft{
    my($self) = @_;
    die if(!exists $self->{'right'});
    die if(!exists $self->{'right'}->{'left'});
    my $root = $self;
    if(defined $self->{'right'}->{'left'}){
	$root = $self->{'right'};
	$self->{'right'} = $self->{'right'}->{'left'};
	$root->{'left'} = $self;
	if(defined $self->{'right'} && defined $self->{'left'}){
	    $self->{'maxend'} = ($self->{'end'}>$self->{'right'}->{'maxend'}) ? $self->{'end'} : $self->{'right'}->{'maxend'};
	    $self->{'maxend'} = ($self->{'maxend'}>$self->{'left'}->{'maxend'}) ? $self->{'maxend'} : $self->{'left'}->{'maxend'};
	    
	    $self->{'minend'} = ($self->{'end'}<$self->{'right'}->{'minend'}) ? $self->{'end'} : $self->{'right'}->{'minend'};
	    $self->{'minend'} = ($self->{'minend'}<$self->{'left'}->{'minend'}) ? $self->{'minend'} : $self->{'left'}->{'minend'};
	}
	elsif(defined $self->{'right'}){
	    $self->{'maxend'} = ($self->{'end'}>$self->{'right'}->{'maxend'}) ? $self->{'end'} : $self->{'right'}->{'maxend'};
	    $self->{'minend'} = ($self->{'end'}<$self->{'right'}->{'minend'}) ? $self->{'end'} : $self->{'right'}->{'minend'};
	}
	elsif(defined $self->{'left'}){ 
	    $self->{'maxend'} = ($self->{'end'}>$self->{'left'}->{'maxend'}) ? $self->{'end'} : $self->{'left'}->{'maxend'};
	    $self->{'minend'} = ($self->{'end'}<$self->{'left'}->{'minend'}) ? $self->{'end'} : $self->{'left'}->{'minend'};
	}
    }
    return $root;
}

1;
