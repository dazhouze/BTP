#!/usr/bin/perl -w
use strict;
my $work_dir = $ARGV[0];
if ($work_dir =~ /(.*)\/$/){
    $work_dir = $1;
}
my @fastqs = split /\s+/, readpipe("ls $work_dir/Fastq | grep fastq");
for my $f (@fastqs){
    $f =~ /phase_([0-1])\.([0-9]+)\.fastq/;
    my $phase = $1;
    my $num = $2;
    print("$phase $num\n");
    my $start = (split/\t/, readpipe("cat $work_dir/Pos/phase_$phase.$num.pos | head -n 1"))[0];
    my $end = (split/\t/, readpipe("cat $work_dir/Pos/phase_$phase.$num.pos | tail -n 1"))[1];
    my $len = $end - $start;
    next if ($len < 1);
    system("canu -pacbio-raw $work_dir/Fastq/phase_$phase.$num.fastq -p frag_$num -d $work_dir/phase_$phase\_$num genomeSize=$len useGrid=false");
}
