#!/usr/bin/perl -w
use strict;
my %list = (
"phase_0.12"=>1,
"phase_0.124"=>1,
"phase_0.134"=>1,
"phase_0.40"=>1,
"phase_0.42"=>1,
"phase_0.57"=>1,
"phase_0.76"=>1,
"phase_0.81"=>1,
"phase_1.39"=>1,
"phase_1.40"=>1,
"phase_1.41"=>1,
"phase_1.42"=>1,
"phase_1.15"=>1,
"phase_1.12"=>1,
"phase_1.133"=>1,
"phase_1.9"=>1,
"phase_1.84"=>1,
"phase_1.57"=>1,
"phase_1.9"=>1,
);
my $work_dir = $ARGV[0];
if ($work_dir =~ /(.*)\/$/){
    $work_dir = $1;
}
my @fastqs = split /\s+/, readpipe("ls $work_dir/Fastq | grep fastq");
for my $f (@fastqs){
    $f =~ /phase_([0-1])\.([0-9]+)\.fastq/;
    my $phase = $1;
    my $num = $2;
    $f =~ /(phase_.*).fastq/;
    my $li = $1;
    if (exists $list{$li}){
        #next if ($phase == 1);
        system("rm -rf canu_assembly/phase_$phase\_$num/*");
        my $start = (split/\t/, readpipe("cat $work_dir/Pos/phase_$phase.$num.pos | head -n 1"))[0];
        my $end = (split/\t/, readpipe("cat $work_dir/Pos/phase_$phase.$num.pos | tail -n 1"))[1];
        my $len = $end - $start;
        next if ($len < 1);
        my $cer = 0.35; # correctedErrorRate
        while ($cer >= 0.25){
            my $contig_num = (split/\s/,readpipe("cat canu_assembly/phase_$phase\_$num/*contigs.fasta | grep \">\" | wc -l "))[0];
            if ($contig_num){
                print("phase_$phase\_$num\t$cer\n") if ($cer != 0.35);
                last;
            }
            else{
                system("rm -rf canu_assembly/phase_$phase\_$num/*");
            }
            system("canu correctedErrorRate=$cer -pacbio-corrected $work_dir/Fastq/phase_$phase.$num.fastq -p frag_$num -d canu_assembly/phase_$phase\_$num genomeSize=$len useGrid=false");
            $cer -= 0.01; # 0.3 0.29 0.28 .. 0.25
        }
    }
}
