#!/usr/bin/perl -w
use strict;

my $mhc_s = 28477797;#start of mhc region of hg19
my $mhc_e = 33448354;#start of mhc region of hg19
$mhc_s -= 1000;#flanking 1k
$mhc_e += 1000;#flanking 1k
my $winSize = 100;#window size (depth)
my $depCut = 10;#depth cutoff value: 2x Canu cutoff value

#BAM file 
my $f ="/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/5YH.best.sorted.bam";
open IN , "samtools view $f | " and print "$f\n";

#resutl file 
my $o = "hardSoftGap.txt";
open OUT , ">$o" and print "Output file: $o\n";

my %depth;#alignment depth
my %depthWin;#hash for depth (win size = 100bp)

while (<IN>){
    next if (/^@/);
    my ($chr,$read_s,$cigar)=(split /\t/)[2,3,5];
    next if ($chr ne "chr6" && $chr ne "6");
    next if ($read_s < 28466797);
    last if ($read_s > $mhc_e);
    #soft clipping at the left most position is removed #$pos is the start of read alignment start
    my $read_e = $read_s;
    while ($cigar =~ /([0-9]+)([A-Z])(.*)/){
        my $num = $1;
        my $typ = $2;
        $cigar = $3;
        if ($typ eq "M" || $typ eq "D"){
            $read_e += $num;
        }
    }
    $read_e -= 1;#
    #depth hash ++;
    for my $p ($read_s .. $read_e){
        $depth{$p}++;
    }
    #hard gap calculation
    for my $w ( (int($read_s/$winSize)+1) .. int($read_e/$winSize)){
        $depthWin{$w}++;
    }
    $depthWin{(int($read_s/$winSize))}+= 1 - ($read_s % $winSize)/$winSize;
    $depthWin{(int($read_e/$winSize)+1)}+= ($read_e % $winSize)/$winSize;
}
close IN ;
print " - Finish read through SAM/BAM file.\n";

#hard gap detect
my $gap_s = 0;
my @hardGap;#2D array for gap region
my $index = 0;
for my $kw (int($mhc_s/$winSize) .. int($mhc_e/$winSize+1)){
    my $value = 0;
    $value = $depthWin{$kw} if (exists $depthWin{$kw});
    if ($value < $depCut && $gap_s == 0){
        $gap_s = $kw*$winSize;
        print OUT "hardS:$gap_s\t"; 
        $hardGap[$index][0] = $gap_s;
    }
    if ($value > $depCut && $gap_s != 0){
        my $dis = $kw*$winSize - $gap_s;
        $gap_s = $kw*$winSize;
        print OUT "hardE:$gap_s\tlen:$dis\n"; 
        $hardGap[$index][1]=$gap_s;
        $gap_s = 0;
        $index++;
    }
}
close OUT;
print " - Finish hard gaps detection.\n";

$index = 0;
for my $kp (($mhc_s+10) .. ($mhc_e-20)){
    next if ($kp > $hardGap[$index][0] && $kp > $hardGap[$index][1]);
    if ($kp>=$hardGap[$index][1]){
        $index++;
        print "go though No.$index hard gap\n";
    }
    #soft gap detect
    #&siginificant == 1: sinificant; else si == 0: not siginificant
    if( &siginificant($kp) ){
        print OUT "soft:$kp\n";
        print "soft:$kp\n";
    }
}
print " - Finish soft gaps detection.\n";

sub siginificant {
    my $pos = $_[0];
    my $n = 10; #sample size
    my @diff;
    for my $kp (($pos-$n+1) .. $pos){
        my $before = 0;
        my $after = 0;
        $before = $depth{$kp} if (exists $depth{$kp});
        $after = $depth{($kp+$n)} if (exists $depth{($kp+$n)});
        push @diff,($before - $after);
    }
    my  $mean_diff = 0;#mean of the diff number list
    for my $val (@diff){
        $mean_diff += $val;
    }
    $mean_diff /= $n;
    my $var_diff = 0;#The variance is a numerical measure of how the data values is dispersed around the mean.
    for my $val (@diff){
        $var_diff +=  ($val-$mean_diff)*($val-$mean_diff);
    }
    $var_diff /= $n;
    my $t_stat_diff = 0;#t score
    if($var_diff){#avoid /0
        $t_stat_diff=abs($mean_diff)/sqrt($var_diff/$n);
    }
    #print "$pos\t$depth{$pos}\t$var_diff\t$t_stat_diff\n";
    if($t_stat_diff > 30){#arbitory cutoff value of t-score
        return 1;
    }
    return 0;
}
