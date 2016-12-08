#!/usr/bin/perl -w
use strict;

my $mhc_s = 28477797;#start of mhc region of hg19
my $mhc_e = 33448354;#start of mhc region of hg19
$mhc_s -= 1000;#flanking 1k
$mhc_e += 1000;#flanking 1k

#BAM file 
my $f ="/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/5YH.best.sorted.bam";
open IN , "samtools view $f | " and print "$f\n";

#resutl file 
my $o = "mhc_flank.depth.txt";
open OUT , ">$o" and print "Output file: $o\n";
print OUT "Coordinate\tDepth\n";

my %depth;#alignment depth

=begin
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
}
close IN ;
=cut

for my $kp (1..50){
    #$depth{$kp} = int(17+rand(3));
    my $p = int 17+rand(3);
    $depth{$kp} = 10*$p;
}
for my $kp (51..200){
    #$depth{$kp} = int(36 + rand(4));
    $depth{$kp} = 10*int(10*int(36 + rand(4))/10);
}

for my $kp ((10) .. (100-10)){
    #print "test: $depth{$kp}\n";
    my $si = &siginificant(\%depth, $kp);
}

=begin
for my $kp (($mhc_s+10) .. ($mhc_e-10)){
    my $value = 0;
    if (exists $depth{$kp}){
        #print OUT "$kp\t$depth{$kp}\n";
        $value = $depth{$kp};
    }   
}
close OUT;
=cut

sub siginificant {
    my ($X, $Y) = @_;
    my %hash = %$X;#depth hash
    my $pos = $Y;
    my $n = 10; #sample size
    my @diff;
    for my $kp (($pos-$n+1) .. $pos){
        push @diff,($hash{$kp}-$hash{($kp+$n)});
    }
    my  $mean_diff = 0;
    for my $val (@diff){
        $mean_diff += $val;
    }
    $mean_diff /= $n;
    my $var_diff = 0;
    for my $val (@diff){
        $var_diff +=  ($val-$mean_diff)*($val-$mean_diff);
    }
    $var_diff /= $n;
    my $t_stat_diff=abs($mean_diff)/sqrt($var_diff/$n);
    #print "$pos\t$depth{$pos}\t$var_diff\t$t_stat_diff\n";
    if($t_stat_diff > 30){
        return 1;
    }
    return 0;
}

