#!/usr/bin/perl -w
use strict;

my $mhc_s = 28477797;#start of mhc region
my $mhc_e = 33448354;#start of mhc region
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

for my $kp ($mhc_s .. $mhc_e){
   if (exists $depth{$kp}){
       print OUT "$kp\t$depth{$kp}\n";
   }
}
close OUT;
