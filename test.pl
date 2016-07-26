#!/usr/bin/perl -w
use strict;
my %hash;#seq depth hash
my $depthCut = 5;

            for my $i (1 .. 50){
                $hash{$i}+=6;
            }
            for my $i (61 .. 100){
                $hash{$i}+=6;
            }
print "finish reading BAM file\n";

my $start;#start coordinary of a gap
my $end;#end coordinary of a gap
my $preGap = 0;#previous gap pos
my $preVal = -1;#previous pos hash value
my $ord = 0;#the num of gaps
for my $i (1 .. 100){
    $hash{$i} = 0 unless (exists $hash{$i});
    print "$i no cover\n" if ($hash{$i} == 0);
    #if($hash{$i} <= $depthCut){#this is some pos called "gap"
    #    if($ord == 0){#the fisrt gap
    #        $ord++;
    #        $start = $i;
    #        $preGap = $i;
    #    }
    #    if ($i == 33448354){#last pos of MHC region hg19
    #        $ord++;
    #        print"gap:$ord\tchr6:$start-$preGap\n";
    #    }
    #    #elsif($i != ($preGap+1)){
    #    elsif($preVal != 0){
    #        $ord++;
    #        print"gap $ord:\tchr6:$start-$preGap\n";
    #        $start =  $i;
    #    }
    #    $preGap = $i;
    #}
    if ($preVal > $depthCut && $hash{$i} <= $depthCut){#gap start
        $ord++;
        print "find gap start\n";
        $start = $i;
    }
    elsif ($preVal <= $depthCut && $hash{$i} > $depthCut){#gap end
        #$ord++;
        $end = $i - 1;
        print"gap $ord:\tchr6:$start-$end\n";
        #$start = $i;
    }
    $preVal = $hash{$i};
}
