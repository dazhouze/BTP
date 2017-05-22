#!/usr/bin/perl -w
use strict;
my $work_dir = $ARGV[0];
if ($work_dir =~ /(.*)\//){
    $work_dir = $1;
}
my $ref_genome = $ARGV[1];
print("$work_dir $ref_genome\n");
my @fastqs = split /\s+/, readpipe("ls $work_dir | grep fastq");
for my $f (@fastqs){
    print("$f\n");
    $f =~ /phase_([0-1])\.([0-9]+)\.fastq/;
    my $phase = $1;
    my $num = $2;
    print("$phase $num\n");
    system("bwa mem -x pacbio $ref_genome $work_dir/$f > $work_dir/phase_$phase.$num.sam");
    system("perl example/SAMbestHit.pl $work_dir/phase_$phase.$num.sam $work_dir/phase_$phase.$num.best.sam");
    system("samtools view -Sb $work_dir/phase_$phase.$num.best.sam > $work_dir/phase_$phase.$num.best.bam");
}
