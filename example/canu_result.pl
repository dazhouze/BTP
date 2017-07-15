#!/usr/bin/perl -w
use strict;
my $work_dir = $ARGV[0];
if ($work_dir =~ /(.*)\/$/){
    $work_dir = $1;
}
my @sub_dir = split /\n/, readpipe("ls $work_dir | grep phase");
for my $f (@sub_dir){
    $f =~ /phase_([0-1])_([0-9]+)/;
    my $phase = $1;
    my $num = $2;

    my $file = readpipe("ls $work_dir/phase_$phase\_$num/");
    if ($file=~/contigs\.fasta/){
        my $contig_num = readpipe("cat $work_dir/phase_$phase\_$num/*contigs.fasta | grep '>' | wc -l");
        chomp($contig_num);
        $contig_num = int($contig_num);
        print("#multi_contig: phase_$phase.$num contig:$contig_num\n") if $contig_num!= 1;
    }else{
        print("#0_contig: phase_$phase.$num\n");
        #system("cat test/canu_assembly_script.sh | grep 'phase_$phase.$num '");
        next;
    }

}
