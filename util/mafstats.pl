#!/usr/bin/perl

#Reports coverage

#Unique DNA should be sum of blocks blocks with one seq and runs aligned to all gaps

use strict;



my $found=0;
my $currscore;
my $currorient;
my $blockorient;
my @allblocks;
my $block = [];
my $isdup=0;
my $multiplealnblkcount=0;
while(my $line=<STDIN>){
    if($line =~ /^a/){
	($currscore) =~ ($line =~ /score=(\S+)/);
	my($label) = ($line =~ /label=(\S+)/);
	my($isdup) = ($line =~ /dup=/) ? 1 : 0;
	push @allblocks,$block if(scalar(@$block)>0);
	$multiplealnblkcount++ if(scalar(@$block)>1);
	$block=[];
    }
    elsif($line =~ /^s/){
	#my @elts = split(/\s+/,$line);
	#0-score,1-blockorient,2-accession,3-start,4-end
	chomp $line;
	push @$block,[$currscore,$line,$isdup];
    }
}
push @allblocks,$block if(scalar(@$block)>0);

#Number of lcbs with N genomes
my $lcbseqcount = [];
#Frequency of alignment columns with N identical rows
my $numIdentCols = [];
#Freq of columns with no gaps
my $numUngappedCols = [];
#Freq columns with one seq and all gaps
my $numGappedCols = [];

#Number of bps in blocks containing N genomes
my $lcbbpcount = [];
my $lcbbpdistro = [];
my $gapdistro = [];


my $alnbpseqs = {};
my $lcbseqs = {};



my $totalscore=0;
my $numgaps=0;
my $numblocks=scalar(@allblocks);
my $totallen=0;
my $totalseqlen=0;
my $smallestblks=0;
my $smallerblks=0;
my $smallerlen=0;
my $smallestlen=0;
my $nummaf=0;

my $uniqcount=0;
my $dupcount=0;

my %minseq;
my %maxseq;
my %allseqs;

print "Num_blocks:$numblocks\n";
print "Num_multi_blocks:$multiplealnblkcount\n";
my $lcbid=0;

open AFILE,"+>aln.$ARGV[0].dat";
foreach my $block (@allblocks){
    my $issmaller=0;
    my $issmallest=0;
    #Min and max len of seqs in the LCB
    my $minlen=-1;
    my $maxlen=0;
    die if(scalar(@$block) ==0);
    my $nseq = scalar(@$block);
    die if($nseq <= 0);
    $lcbseqcount->[$nseq]++;
    my @alntext;
    my $isdup=0;
    if($nseq>1){
	foreach my $maf (@$block){
	    if($maf->[2]){
		$isdup=1;
	    }
	    my($seq,$beg,$len,$orient,$seqlen,$text) = ($maf->[1] =~ /s\s+(\S+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+)\s+(\S+)/);
	    die if($len<0);
	    $text =~ s/\s+//g;
	    if($text =~ /[^-]/){
		push @alntext,$text;
		if(exists $allseqs{$seq}){
		    die if($seqlen != $allseqs{$seq});
		}
		else{
		    $allseqs{$seq} = $seqlen;
		}
		if($minlen==-1){
		    $minlen=$len;
		}
		else{
		    $minlen = ($len<$minlen) ? $len:$minlen;
		    die if($minlen<0);
		}
		$maxlen = ($len>$maxlen) ? $len:$maxlen;
		my $cgaps = ($text =~ tr/\-/-/);
		die if($cgaps<0);
		$numgaps += $cgaps;
		my($fmin,$fmax);
		
		if($orient eq '-'){
		    $fmin = $seqlen-$beg-$len;
		    $fmax = $seqlen-$beg;
		}
		else{
		    $fmin = $beg;
		    $fmax = $beg+$len;
		}
		die "$maf->[1]" if($fmin < 0 || $fmin > $seqlen);
		die "$maf->[1]" if($fmax < 0 || $fmax > $seqlen);
		$minseq{$seq} = ($minseq{$seq} < $fmin) ? $fmin : $minseq{$seq};
		$maxseq{$seq} = ($maxseq{$seq} > $fmax) ? $fmax : $maxseq{$seq};

		$lcbseqs->{$seq} = [] if(!ref $lcbseqs->{$seq});
		push @{$lcbseqs->{$seq}},[$fmin,$fmax,$orient];

		die "$maf->[1]" if($len<=0);
		$totallen += $len;
		$nummaf++;
		if($len < 100){
		    $issmallest=1;
		    $smallestlen+=$len;
		}
		if($len < 1000){
		    $issmaller=1;
		    $smallerlen+=$len;
		}
	    }
	    else{
		print STDERR "All gap encountered but length $len > 0 $text\n" if($len != 0);
		$nseq--;
	    }
	}
    }
    else{
	my($seq,$beg,$len,$orient,$seqlen,$text) = ($block->[0]->[1] =~ /s\s+(\S+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+)\s+(\S+)/);
	$minlen=$len;
    }
    die if($minlen<0);
    $smallerblks++ if($issmaller);
    $smallestblks++ if($issmallest);

    $lcbbpcount->[$nseq]+=$minlen;

    if($nseq ==1){
	if($isdup){
	    $uniqcount +=$minlen;
	}
	else{
	    $dupcount +=$minlen;
	}
    }
	    
    $lcbbpdistro->[$nseq] = [] if(!ref $lcbbpdistro->[$nseq]);
    push @{$lcbbpdistro->[$nseq]},$minlen;

    print STDERR "LCB: $lcbid maxlen:$maxlen\n";
    $lcbid++;

    if($nseq>1){
	my $alnmatrix = &maf2matrix(\@alntext);
	my($lcbtotalscore,$blklen) = &scorealn($alnmatrix,
					       $numIdentCols,
					       $numUngappedCols,
					       $numGappedCols,
					       $gapdistro,
					       $alnbpseqs);
	die "$lcbtotalscore,$blklen" if($blklen ==0);
	$totalscore+=$lcbtotalscore;
	my $alnlen = scalar(@{$alnmatrix->[0]});
	my $nseq = scalar(@$alnmatrix);
	for(my $k=0;$k<$alnlen;$k++){
	    for(my $i=0;$i<$nseq;$i++){
		print AFILE $alnmatrix->[$i]->[$k];
	    }
	    print AFILE "\n";
	}
    }
}

close AFILE;
my $seqmatrix = &getCovered($lcbseqs);

my $uniqbptotal=0;
my $alignedlentotal=0;
my $doublecovtotal=0;
open MFILE,"+>bps.$ARGV[0].dat";
foreach my $seq (sort {$a cmp $b} keys %allseqs){
    $totalseqlen+=$allseqs{$seq};
    my $alignedlen=0;
    my $doublecov=0;
    my $uniqbp=0;
    for(my $i=0;$i<$allseqs{$seq};$i++){    
	if($seqmatrix->{$seq}->[$i]){
	    $alignedlen++;
	    $alignedlentotal++;
	    if($seqmatrix->{$seq}->[$i]>1){
		$doublecovtotal+=$seqmatrix->{$seq}->[$i]-1;
		$doublecov++;
	    }
	}
	else{
	    die if($seqmatrix->{$seq}->[$i]>0);
	    $uniqbptotal++;
	    $uniqbp++;
	}
	print MFILE "$seq $i $seqmatrix->{$seq}->[$i]\n";
    }
    print "$seq len:$allseqs{$seq} aln_cov:$alignedlen aln_cov_pct:",$alignedlen/$allseqs{$seq}," uniq:$uniqbp doublecov:$doublecov \n";
}
close MFILE;

#Count of bases that are aligned to only gaps
my $uniqaln=0;
for(my $i=0;$i<scalar(@$numGappedCols);$i++){
    $uniqaln+=$numGappedCols->[$i];
}

print "\n";
#Summary #genomes,total len, avg block size
print "max_genomes_aln:",scalar(@$numIdentCols)-1,"\n";
print STDERR "Num ident cols size=",scalar(@$numIdentCols),"!= Numbpdistro=",scalar(@$lcbbpdistro),"\n" if(scalar(@$numIdentCols)!=scalar(@$lcbbpdistro));
print "total_seq_len:",$totalseqlen,"\n";
print "avg_block_len:",$totallen/$nummaf,"\n";
print "num_lcbs:",$nummaf,"\n";
print "double_covered:",$doublecovtotal,"\n";
#Avg/total coverage, #bps aligned
print "aln_cov:",$alignedlentotal," ",$totallen-$doublecovtotal,"\n";
print "aln_cov_pct:",$alignedlentotal/$totalseqlen,"\n";
print "not_cov:",$uniqbptotal,"\n";
print "not_cov_pct:",$uniqbptotal/$totalseqlen,"\n";

#Composition
print "aln_bps:",($totalseqlen-$uniqbptotal-$uniqaln),"\n";
print "aln_pct:",($totalseqlen-$uniqbptotal-$uniqaln)/$totalseqlen,"\n";
print "core_bps:",$numUngappedCols->[scalar(@$numUngappedCols)-1],"\n";
print "core_pct:",$numUngappedCols->[scalar(@$numUngappedCols)-1]/$totalseqlen,"\n";
print "uniq_bps:",$uniqbptotal+$uniqaln,"\n";
print "uniq_pct:",($uniqbptotal+$uniqaln)/$totalseqlen,"\n";

print "\n";
print "MISMATCH between uniqLCB len and calculated len\n" if( $lcbbpcount->[1] != $uniqbptotal);
print "uniq_LCBlen:",$lcbbpcount->[1],"\n";
print "uniq_cov:",$uniqbptotal,"\n";
print "uniq_aln:",$uniqaln,"\n";
print "uniq_dup:",$uniqcount,"\n";
print "dup_bps:",$dupcount,"\n";

print "blklt100bp:",$smallestblks,"\n";
print "blklen:",$smallestlen,"\n";
print "blklt1000bp:",$smallerblks,"\n";
print "blklen:",$smallerlen,"\n";
#Scoring
print "num_gaps:",$numgaps,"\n";
print "score:",$totalscore,"\n";


print "LCB seq count\n";
for(my $i=0;$i<scalar(@$lcbseqcount);$i++){
    print "$i\t";
}
print "\n";
for(my $i=0;$i<scalar(@$lcbseqcount);$i++){
    print $lcbseqcount->[$i],"\t";
}
print "\n";
print "LCB coverage bp count\n";
for(my $i=0;$i<scalar(@$lcbbpcount);$i++){
    print "$i\t";
}
print "\n";
for(my $i=0;$i<scalar(@$lcbbpcount);$i++){
    print $lcbbpcount->[$i],"\t";
}
print "\n";
print "Ident.Freq of identical alignment columns\n";
for(my $i=0;$i<scalar(@$numIdentCols);$i++){
    print "$i\t";
}
print "\n";
for(my $i=0;$i<scalar(@$numIdentCols);$i++){
    print "$numIdentCols->[$i]\t";
}
print "\n";
print "NoGaps.Freq of alignment columns with no gaps\n";
for(my $i=0;$i<scalar(@$numUngappedCols);$i++){
    print "$i\t";
}
print "\n";
for(my $i=0;$i<scalar(@$numUngappedCols);$i++){
    print "$numUngappedCols->[$i]\t";
}
print "\n";
print "AllGaps.Freq of alignment cols with one seq and all gaps\n";
for(my $i=0;$i<scalar(@$numGappedCols);$i++){
    print "$i\t";
}
print "\n";
for(my $i=0;$i<scalar(@$numGappedCols);$i++){
    print "$numGappedCols->[$i]\t";
}
print "\n";

print "LCBs:";
my @lcblens;
for(my $i=2;$i<@$lcbbpdistro;++$i){
    push @lcblens,@{$lcbbpdistro->[$i]} if(ref $lcbbpdistro->[$i]);
}
print join(',',sort {$a <=> $b} @lcblens);
print "\n";
print "LCBs core:";
print join(',',sort {$a <=> $b} @{$lcbbpdistro->[scalar(@$lcbbpdistro)-1]});
print "\n";
print "Gaps:";
my @gaplens;
foreach my $seq (@$gapdistro){
    push @gaplens,@$seq if(ref $seq);
}
print "\n";
print join(',',sort {$a <=> $b} @gaplens);
print "\n";

print STDERR "Writing .dat files for R\n";

#Data for R
open LFILE,"+>lcbs.$ARGV[0].dat";
print LFILE join("\n",sort {$a <=> $b} @lcblens);
close LFILE;

open CFILE,"+>corelcbs.$ARGV[0].dat";
print CFILE join("\n",sort {$a <=> $b} @{$lcbbpdistro->[scalar(@$lcbbpdistro)-1]});
close CFILE;

open GFILE,"+>gaps.$ARGV[0].dat";
print GFILE join("\n",sort {$a <=> $b} @gaplens);
close GFILE;

open RFILE,"+>mafstats.$ARGV[0].r";
print RFILE "lcbs <- read.csv(file=\"lcbs.$ARGV[0].dat\");\n";
print RFILE "corelcbs <- read.csv(file=\"corelcbs.$ARGV[0].dat\");\n";
print RFILE "gaps <- read.csv(file=\"gaps.$ARGV[0].dat\");\n";
print RFILE "hist(lcbs\$X1, col=\"green\", main=\"LCBs\", xlab=\"LCB length (bp)\");\n";
print RFILE "dev.print(device=postscript, \"lcbs.$ARGV[0].eps\", onefile=FALSE, horizontal=FALSE);\n";
print RFILE "hist(corelcbs\$X1, col=\"blue\", main=\"Core LCBs\", xlab=\"LCB length (bp)\");\n";
print RFILE "dev.print(device=postscript, \"corelcbs.$ARGV[0].eps\", onefile=FALSE, horizontal=FALSE);\n";
print RFILE "hist(gaps\$X1, col=\"red\", main=\"Gaps\", xlab=\"Gap length (bp)\");\n";
print RFILE "dev.print(device=postscript, \"gaps.$ARGV[0].eps\", onefile=FALSE, horizontal=FALSE);\n";
close RFILE;

sub scorealn{
    my($matrix,$numIdentCols,$numUngappedCols,$numGappedCols,$gapaln,$alnbpseqs) = @_;
    my $gapext = -1;
    my $gapopen = -2;
    my $gapopeni=0;
    my $gapopenj=0;
    my $gapexcount=0;
    my $gapcount=0;
    my $totalscore=0;
    my $alnlen = 0;
    my $nseq = scalar(@$matrix);
    #print "Scoring $nseq\n";

    #Loop over each sequence/row
    for(my $i=0;$i<$nseq;$i++){
	if($alnlen!=0){
	    die if($alnlen != scalar(@{$matrix->[$i]}));
	}
	else{
	    $alnlen = scalar(@{$matrix->[$i]});
	}
	for(my $j=$i+1;$j<$nseq;$j++){
	    
	    #print "$i $alnlen\n";
	    #Loop over each column
	    for(my $k=0;$k<$alnlen;$k++){
		if($matrix->[$i]->[$k] ne '-'){
		    if($matrix->[$j]->[$k] ne '-'){
			$gapopeni=0;
			$gapopenj=0;
			$totalscore+=1;#$scorematrix[$matrix[$i][$k]][$matrix[$j][$k]];
		    }
		    else{
			if($gapopenj){
			    $gapexcount++;
			    $totalscore+=$gapext;
			}
			else{
			    $gapopenj=1;
			    $gapcount++;
			    $totalscore+=$gapopen;
			}
		    }
		}
		else{
		    if($matrix->[$j]->[$k] ne '-'){
			if($gapopeni){
			    $gapexcount++;
			    $totalscore+=$gapext;
			}
			else{
			    $gapopeni=1;
			    $gapcount++;
			    $totalscore+=$gapopen;
			}
		    }
		}
	    }
	}
    }
    #Get number of identical columns, allowing for gaps but not mismatches
    # S1 TTTTTTAAATTT
    # S2 TT---TAAAA-A
    # S3 TTTTTT--ATTT
    #    332223223020
    # $numIdentCols[0]=2 //at least one mismatch
    # $numIdentCols[2]=6
    # $numIdentCols[3]=3
    my $c; #bp
    my @uniqruns;
    my $uniqrow;
    my $runopen;
    my $startrun=-1;
    my $runpos=-1;
    for(my $k=0;$k<$alnlen;$k++){
	my $numIdents=0;
	my $mismatch;
       	for(my $j=0;$j<$nseq;$j++){
	    if($matrix->[$j]->[$k] ne '-'){
		if($numIdents==0){
		    $c = lc($matrix->[$j]->[$k]);
		    $numIdents++;
		    $uniqrow=$j;
		}
		else{
		    if(lc($matrix->[$j]->[$k]) eq $c){
			$numIdents++;
		    }
		    else{
			$numIdents=0;
			last;
		    }
		}
	    }
	    else{
		$mismatch=1;
	    }
	}
	if($numIdents==1){
	    if($runopen eq $uniqrow){
		$runpos=$k;
	    }
	    else{
		push @uniqruns,[$runopen,$startrun,$runpos] if($runopen ne "");
		$runopen=$uniqrow;
		$startrun=$k;
		$runpos=$k;
	    }
	}
	else{
	    push @uniqruns,[$runopen,$startrun,$runpos] if($runopen ne "");
	    $runopen = "";
	    $startrun=-1;
	    $runpos=-1;
	}
	$numIdentCols->[$numIdents]++;
	$mismatch =1 if($numIdents<$nseq);
	#push @$mismatches,$k if($mismatch);
    }
    #Get number of ungapped columns
    # S1 TTTTTTAAATTT
    # S2 TT---TAAAA-A
    # S3 TTTTTT--ATTT
    #    330003003303
    # $numUngapped[0]=6
    # $numUngapped[3]=6
    my $c; #bp
    for(my $k=0;$k<$alnlen;$k++){
	my $numUngaps=0;
	for(my $j=0;$j<$nseq;$j++){
	    if($matrix->[$j]->[$k] ne '-'){
		$numUngaps++;
	    }
	    else{
		$numUngaps=0;
		last;
	    }
	}
	$numUngappedCols->[$numUngaps]++;
    }
    #Get number columns with one sequence and all gaps
    for(my $k=0;$k<$alnlen;$k++){
	my $numGaps=0;
	for(my $j=0;$j<$nseq;$j++){
	    if($matrix->[$j]->[$k] eq '-'){
		$numGaps++;
	    }
	}
	if($numGaps==$nseq){
	    for(my $j=0;$j<$nseq;$j++){
		for(my $k=0;$k<$alnlen;$k++){
		    print STDERR "$matrix->[$j]->[$k]";
		}
		print STDERR "\n";
	      }  
	    print STDERR "Column $k has all gaps\n";
	}
	if($numGaps>0 && $numGaps==$nseq-1){
	    $numGappedCols->[$numGaps]++;
	}
    }

    #
    # Save lengths of all runs of gaps
    # Eg.
    # S1 TTTTTTAAATTT
    # S2 TT---TAAAT-T
    # S3 TTTTTT--ATTT
    # 
    #Results
    #$gapaln[1]=[3,1]
    #$gapaln[2]=[2]
    $gapopen=0;    
    for(my $j=0;$j<$nseq;$j++){
	$gapopen=0;
	for(my $k=0;$k<$alnlen;$k++){
	    if($matrix->[$j]->[$k] eq '-'){
		$gapopen++;
	    }
	    else{
		if($gapopen){ #end of a run of gaps
		    $gapaln->[$j] = [] if(!ref $gapaln->[$j]);
		    push @{$gapaln->[$j]},$gapopen;
		    #if($gapopen>1000){
		    #print STDERR "Long gap $gapopen in seq $j\n";
		    #}
		}
		#start of a run of gaps
		$gapopen=0;
	    }
	}
    }
    if($gapopen){
	$gapaln->[$nseq-1] = [] if(!ref $gapaln->[$nseq-1]);
	push @{$gapaln->[$nseq-1]},$gapopen;
    }
    return ($totalscore,$alnlen);
}

sub maf2matrix{
    my($mafs) = @_;
    my $matrix = [];
    my $i=0;
    print STDERR " with ",scalar(@$mafs)," seqs\n";
    foreach my $m (@$mafs){
	my @row = split(//,$m);
	$matrix->[$i++] = \@row;
    }
    return $matrix;
}



sub getCovered{
    my($blocksbyseq) = @_;
    my $seqmatrix = {};
    
    foreach my $seq (sort {$a cmp $b} keys %$blocksbyseq){
	foreach my $b (@{$blocksbyseq->{$seq}}){
	    for(my $j=$b->[0];$j<$b->[1];$j++){
		if($seqmatrix->{$seq}->[$j]>0){
		    print STDERR " $seq $j doublecov $seqmatrix->{$seq}->[$j] $b->[0] $b->[1]\n";
		}
		$seqmatrix->{$seq}->[$j]++;
	    }
	}
    }
    return $seqmatrix;
}    
	    

