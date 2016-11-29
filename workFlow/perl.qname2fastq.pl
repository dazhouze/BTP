#!/usr/bin/perl -w
use strict;
my $fq="/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS.fastq";
for my $f (@ARGV){
    my %hash;
    open IN , "$f";
    while (<IN>){
        chomp;
        my $qname=(split /\s+/)[0];
        $hash{$qname}++;
    }
    close IN;
    open IN2,"$fq" or die"$fq $?\n";
    open OUT , ">$f.fastq";
    while (<IN2>){
        chomp;
        my $one=$_;
        my $two=<IN2>;
        my $thr=<IN2>;
        my $fou=<IN2>;
        $one=~/@(.*)$/;
        my $q=$1;
        if (exists $hash{$q}){
            print OUT "$one\n$two$thr$fou";
        }
    }
    close IN2;
    close OUT;
}
