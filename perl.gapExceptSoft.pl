#!/usr/bin/perl -w
use strict;
my %hash;#seq depth hash

open OUT , ">gapExcepthSoftClipping.txt";
my $time = localtime();
print "Start in $time\n";
open IN , "/home/zhouze/team/YH05_PacBio/YH06_YH07/2_bestHit/YH14567.best.sam";
while(<IN>){
    next if (/^@/);
    chomp;
    my $con = $_;
    next unless (&ifMHC($con));#MHC region overlap
    &seqCover($con, \%hash);
}
close IN;
$time = localtime();
print "Finish reading BAM file in $time\n";
for my $cutoffValue (1 .. 10){
    &cufoff($cutoffValue);
}
$time = localtime();
print "Finish traverse all hash in $time\n";
close OUT;

########## ########## Functions ########## ##########
sub cufoff {
    my $depthCut =  $_[0];
    my @start;#start coordinary of a gap
    my $end;#end coordinary of a gap
    my $preGap = 0;#previous gap pos
    my $preVal = $depthCut + 1;#previous pos hash value
    my $ord = 0;#the num of gaps
    for my $i (28477797 .. 33448354){
        $hash{$i} = 0 unless (exists $hash{$i});
        #print "$i no cover\n" if ($hash{$i} == 0);
        if ($preVal > $depthCut && $hash{$i} <= $depthCut){#gap start
            $ord++;
            #print "find gap start\n";
            push @start, $i;
        }
        elsif ($preVal <= $depthCut && $hash{$i} > $depthCut){#gap end
            $end = $i - 1;
            my $st = pop @start;
            print OUT "cutoff: $depthCut gap $ord:\tchr6:$st-$end\n";
        }
        $preVal = $hash{$i};
    }
    if($#start > 0){#the end of MHC region is a gap position
        my $st = pop @start;
        $end = 33448354 - 1;
        print OUT "cutoff: $depthCut gap $ord:\tchr6:$st-$end\n";
    }
}

sub ifMHC {
    my $con = shift;
    my @col = split /\t/, $con;
    if ($col[2] eq "chr6"){
        unless ($col[3]>33448354 ){
            unless (($col[3]+length($col[9]))<28477797){
                return 1;
            }
        }
    }
    return 0;
}

sub seqCover {
    my $con = $_[0];
    my $X = $_[1];#hash

    my ($qname,$chr,$pos,$cigar,$seq,$qual)=(split /\t/,$con)[0,2,3,5,9,10];
    my $now_ci=$cigar;
    my $now_pos=$pos;#position in reference genome
    
    #use cigar field to figure out base N.O, base, type.
    $now_pos=$pos-1;
    while($now_ci=~ /^([0-9]+)([A-Z])(.*)/){
        my $num=$1;
        my $type=$2;
        $now_ci=$3;
        if ($type eq "M"){
            for my $i (1 .. $num){
                $now_pos++;
                $$X{$now_pos}++;
            }
        }elsif ($type eq "I"){
        }elsif ($type eq "D"){
                $now_pos+=$num;
        }elsif ($type eq "S"){
        }elsif ($type eq "H"){
        }else{
            print "Unknow Cigar character.\n";
            die;
        }
    }
    return 0;
}
