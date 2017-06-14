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
    next if ($phase == 0);
    my $start = (split/\t/, readpipe("cat $work_dir/Pos/phase_$phase.$num.pos | head -n 1"))[0];
    my $end = (split/\t/, readpipe("cat $work_dir/Pos/phase_$phase.$num.pos | tail -n 1"))[1];
    my $len = $end - $start;
    next if ($len < 1);
    my $cer = 0.3; # correctedErrorRate
    while ($cer >= 0.04){
        my $contig_num = (split/\s/,readpipe("cat canu_assembly/phase_$phase\_$num/*contigs.fasta | grep \">\" | wc -l "))[0];
        if ($contig_num){
            print("phase_$phase\_$num\t$cer\n");
            last;
        }
        else{
            system("rm -rf canu_assembly/phase_$phase\_$num/*");
        }
        system("canu correctedErrorRate=$cer -pacbio-corrected $work_dir/Fastq/phase_$phase.$num.fastq -p frag_$num -d canu_assembly/phase_$phase\_$num genomeSize=$len useGrid=false");
        $cer -= 0.005; # 0.3 0.29 0.28 .. 0.25
    }
}
