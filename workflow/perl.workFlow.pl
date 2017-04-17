#!/usr/bin/perl -w
use strict;
use Getopt::Std;

print"\nNotice:\n\tGenes coordinates are based on hg19!HLA-DPA1 and HLA-DPB1 merged into one directory\n";
my %opts;
getopts('s:f:q:o:', \%opts);
$opts{o}="./" unless ($opts{o});
print"BestHit SAMFile:$opts{s}\nReferenceGenome:$opts{f}\nFastqFile:$opts{q}\nOutputDirectory:$opts{o}\n";
die "Usage:\n\t-f for reference fasta.\n\t-o for output directory\n\t-q for fastq files.\n\t-s for SAM files.\n" unless ($opts{s} && $opts{f} && $opts{q} && $opts{o});
my $sam=$opts{s};
my $fq=$opts{q};
my $outDir=$opts{o};
system "mkdir -p $outDir";
system "mkdir -p $outDir/SAM";
system "mkdir -p $outDir/BAM";
system "mkdir -p $outDir/Raw";
system "mkdir -p $outDir/Phase";
system "mkdir -p $outDir/Result";
system "mkdir -p $outDir/Temp";

my %gene;
$gene{"A"}=([29910247,29913661]);
$gene{"B"}=([31321649,31324989]);
$gene{"C"}=([31236526,31239913]);
$gene{"DQA1"}=(["32605169","32612152"]);
$gene{"DQB1"}=(["32627241","32634466"]);
$gene{"DPAB1"}=(["33032346","33057473"]);
$gene{"DRB1"}=(["32546546","32557613"]);

my @bam=qw"A.0.bam A.1.bam B.0.bam B.1.bam C.0.bam C.1.bam DPAB1.0.bam DPAB1.1.bam  DQA1.0.bam DQA1.1.bam DQB1.0.bam DQB1.1.bam  DRB1.0.bam  DRB1.1.bam ";

################### fetch sam ccs read from fastq (mapq>=20 len>= 1000)
for my $kgene (sort keys %gene){
    my $start=$gene{$kgene}[0];
    my $end=$gene{$kgene}[1];
    my %hash;
    open IN , "$sam" or die "$sam $?\n";
    my $read=0;
    my $base=0;
    while (<IN>){
        chomp;
        next if (/^@/);
        my ($qname,$chr,$pos,$mapq,$seq,$qual)=(split /\t/)[0,2,3,4,9,10];
        my $len=length($seq);
        if ($chr eq "chr6" && $mapq>=20 && $len>= 1000){
            my $end0=$pos + length($seq);
            if ($pos>$end || $end0<$start){
                next;
            }elsif ($pos<$start && $end0>$start){
                $read++;
                $base+=length $seq;
                $hash{$qname}++ ;
            }
            elsif ($pos<$end && $end0>$end){
                $read++;
                $base+=length $seq;
                $hash{$qname}++ ;
            }
            elsif ($pos>$start && $end0<$end){
                $read++;
                $base+=length $seq;
                $hash{$qname}++ ;
            }
        }
    }
    close IN;

    open IN2,"$fq" or die"$fq $?\n";
    open OUT , ">$outDir/Raw/$kgene.fastq";
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
    print"Finish extract HLA-$kgene\n";
}

##################### best hit
sub bestHit{
    open IN , "$_[0]";
    open OUT , ">$_[1]";
    my %hash;#store read align
    while (<IN>){
        print OUT "$_" if (/^@/);
        next if (/^@/);
        chomp;
        $_=~/AS:i:(.*?)\s/;
        my $as=$1;#AS field
        my ($qname)=(split /\t/)[0];#QNAME to identify read
        $hash{$qname}{$as}="$_";
    }
    close IN;
    #find best hit align
    for my $kq (keys %hash){
        my @AS=reverse sort {$a<=>$b} keys %{$hash{$kq}};
        print OUT "$hash{$kq}{$AS[0]}\n";
    }
    close OUT;
}

##################### bwa mem for phasing
&alignment;
sub alignment{
    system"bwa mem $opts{f} $outDir/Raw/A.fastq    > $outDir/SAM/A.sam";
    system"bwa mem $opts{f} $outDir/Raw/B.fastq    > $outDir/SAM/B.sam";
    system"bwa mem $opts{f} $outDir/Raw/C.fastq    > $outDir/SAM/C.sam";
    system"bwa mem $opts{f} $outDir/Raw/DPAB1.fastq > $outDir/SAM/DPAB1.sam";
    system"bwa mem $opts{f} $outDir/Raw/DQA1.fastq > $outDir/SAM/DQA1.sam";
    system"bwa mem $opts{f} $outDir/Raw/DQB1.fastq > $outDir/SAM/DQB1.sam";
    system"bwa mem $opts{f} $outDir/Raw/DRB1.fastq > $outDir/SAM/DRB1.sam";
}

&bestHit   ("$outDir/SAM/A.sam","$outDir/SAM/A.bestHit.sam");
&bestHit   ("$outDir/SAM/B.sam","$outDir/SAM/B.bestHit.sam");
&bestHit   ("$outDir/SAM/C.sam","$outDir/SAM/C.bestHit.sam");
&bestHit("$outDir/SAM/DPAB1.sam","$outDir/SAM/DPAB1.bestHit.sam");
&bestHit("$outDir/SAM/DQA1.sam","$outDir/SAM/DQA1.bestHit.sam");
&bestHit("$outDir/SAM/DQB1.sam","$outDir/SAM/DQB1.bestHit.sam");
&bestHit("$outDir/SAM/DRB1.sam","$outDir/SAM/DRB1.bestHit.sam");

&sortIndex;
sub sortIndex{
    system"samtools view -Sb $outDir/SAM/A.bestHit.sam    | samtools sort > $outDir/BAM/A.sort.bam";    
    system"samtools view -Sb $outDir/SAM/B.bestHit.sam    | samtools sort > $outDir/BAM/B.sort.bam";    
    system"samtools view -Sb $outDir/SAM/C.bestHit.sam    | samtools sort > $outDir/BAM/C.sort.bam";    
    system"samtools view -Sb $outDir/SAM/DPAB1.bestHit.sam | samtools sort > $outDir/BAM/DPAB1.sort.bam"; 
    system"samtools view -Sb $outDir/SAM/DQA1.bestHit.sam | samtools sort > $outDir/BAM/DQA1.sort.bam"; 
    system"samtools view -Sb $outDir/SAM/DQB1.bestHit.sam | samtools sort > $outDir/BAM/DQB1.sort.bam"; 
    system"samtools view -Sb $outDir/SAM/DRB1.bestHit.sam | samtools sort > $outDir/BAM/DRB1.sort.bam"; 

    system"samtools index $outDir/BAM/A.sort.bam";
    system"samtools index $outDir/BAM/B.sort.bam";
    system"samtools index $outDir/BAM/C.sort.bam";
    system"samtools index $outDir/BAM/DPAB1.sort.bam";
    system"samtools index $outDir/BAM/DQA1.sort.bam";
    system"samtools index $outDir/BAM/DQB1.sort.bam";
    system"samtools index $outDir/BAM/DRB1.sort.bam";
}

#&phase;
sub phase{
    system"samtools phase -b $outDir/Temp/A    $outDir/BAM/A.sort.bam";   
    system"samtools phase -b $outDir/Temp/B    $outDir/BAM/B.sort.bam";   
    system"samtools phase -b $outDir/Temp/C    $outDir/BAM/C.sort.bam";   
    system"samtools phase -b $outDir/Temp/DPAB1 $outDir/BAM/DPAB1.sort.bam";
    system"samtools phase -b $outDir/Temp/DQA1 $outDir/BAM/DQA1.sort.bam";
    system"samtools phase -b $outDir/Temp/DQB1 $outDir/BAM/DQB1.sort.bam";
    system"samtools phase -b $outDir/Temp/DRB1 $outDir/BAM/DRB1.sort.bam";
}

##################### phasing check
my @zeroFile=split /\n/,readpipe"ls $outDir/Temp/ | gerp 0.bam";
my @oneFile=split /\n/,readpipe"ls $outDir/Temp/ | gerp 1.bam";

##################### BAM to fastq
for my $bam (@bam){
    my %hash;
    open IN , "samtools view $outDir/Phase/$bam |";
    $bam=~/(.*)bam/;
    my $out="$1"."fastq";
    open OUT , ">$outDir/Phase/$out";
    #print"bam $bam\tfastq $out\n";
    my $read=0;
    my $base=0;
    while (<IN>){
        chomp;
        my ($qname,$chr,$pos,$seq,$qual)=(split /\t/)[0,2,3,9,10];
        if ($chr eq "chr6"){
                $read++;
                $base+=length $seq;
                $hash{$qname}++;
        }
    }
    close IN;

    open IN2,"$fq";
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
######################### Canu assembly
#&assembly;
sub assembly{
    system"canu  useGrid=false -p A -d $outDir/A   genomesize=5000 -pacbio-raw  $outDir/Phase/A.fastq";
    system"canu  useGrid=false -p A -d $outDir/A_0 genomesize=5000 -pacbio-raw  $outDir/Phase/A.0.fastq";
    system"canu  useGrid=false -p A -d $outDir/A_1 genomesize=5000 -pacbio-raw  $outDir/Phase/A.1.fastq";
    system"canu  useGrid=false -p B -d $outDir/B   genomesize=5000 -pacbio-raw  $outDir/Phase/B.fastq";
    system"canu  useGrid=false -p B -d $outDir/B_0 genomesize=5000 -pacbio-raw  $outDir/Phase/B.0.fastq";
    system"canu  useGrid=false -p B -d $outDir/B_1 genomesize=5000 -pacbio-raw  $outDir/Phase/B.1.fastq";
    system"canu  useGrid=false -p C -d $outDir/C   genomesize=5000 -pacbio-raw  $outDir/Phase/C.fastq";
    system"canu  useGrid=false -p C -d $outDir/C_0 genomesize=5000 -pacbio-raw  $outDir/Phase/C.0.fastq";
    system"canu  useGrid=false -p C -d $outDir/C_1 genomesize=5000 -pacbio-raw  $outDir/Phase/C.1.fastq";
    system"canu  useGrid=false -p DPAB1 -d $outDir/DPAB1   genomesize=7000 -pacbio-raw    $outDir/Phase/DPAB1.fastq";
    system"canu  useGrid=false -p DPAB1 -d $outDir/DPAB1_0 genomesize=7000 -pacbio-raw    $outDir/Phase/DPAB1.0.fastq";
    system"canu  useGrid=false -p DPAB1 -d $outDir/DPAB1_1 genomesize=7000 -pacbio-raw    $outDir/Phase/DPAB1.1.fastq";
    system"canu  useGrid=false -p DQA1 -d $outDir/DQA1   genomesize=7000 -pacbio-raw    $outDir/Phase/DQA1.fastq";
    system"canu  useGrid=false -p DQA1 -d $outDir/DQA1_0 genomesize=7000 -pacbio-raw    $outDir/Phase/DQA1.0.fastq";
    system"canu  useGrid=false -p DQA1 -d $outDir/DQA1_1 genomesize=7000 -pacbio-raw    $outDir/Phase/DQA1.1.fastq";
    system"canu  useGrid=false -p DQB1 -d $outDir/DQB1   genomesize=6000 -pacbio-raw    $outDir/Phase/DQB1.fastq";
    system"canu  useGrid=false -p DQB1 -d $outDir/DQB1_0 genomesize=6000 -pacbio-raw    $outDir/Phase/DQB1.0.fastq";
    system"canu  useGrid=false -p DQB1 -d $outDir/DQB1_1 genomesize=6000 -pacbio-raw    $outDir/Phase/DQB1.1.fastq";
    system"canu  useGrid=false -p DRB1 -d $outDir/DRB1   genomesize=5000 -pacbio-raw    $outDir/Phase/DRB1.fastq";
    system"canu  useGrid=false -p DRB1 -d $outDir/DRB1_0 genomesize=5000 -pacbio-raw    $outDir/Phase/DRB1.0.fastq";
    system"canu  useGrid=false -p DRB1 -d $outDir/DRB1_1 genomesize=5000 -pacbio-raw    $outDir/Phase/DRB1.1.fastq";
}
