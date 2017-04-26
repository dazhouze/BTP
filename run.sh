#!/bin/bash
dir=./test
rm -rf $dir/
# low long region 91k
#python phasing.py -s 28497412   -e 28588898
# A gene
python phasing.py -s 29910247   -e 29913661 -o $dir 
# 7 k
#python phasing.py -s 29053005 -e 29060759 -o $dir
cp $dir/log/* ~/D*
exit


python fetch.py $dir/phase_0.txt $dir/phase_1.txt
# BWA Alignment Command
align="bwa mem -x pacbio"
# Homo Sapiens Reference Renome hg19
refe=/ifs1/ST_IM/USER/database/Resequence/database/Ref/hg19/ucsc.hg19.fa
# SAM File Best Hit Perl Script
best="/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/workflow/step_2_CCS2BAM/SAMbestHit.pl" 
# SAM To BAM Format Conversion Command
trans="samtools view -Sb"
$align $refe $dir/phase_0.fq > $dir/phase_0.sam
$align $refe $dir/phase_1.fq > $dir/phase_1.sam
perl $best $dir/phase_0.sam $dir/phase_0.best.sam
perl $best $dir/phase_1.sam $dir/phase_1.best.sam
$trans $dir/phase_0.best.sam > $dir/phase_0.best.bam
$trans $dir/phase_1.best.sam > $dir/phase_1.best.bam
samtools sort -o $dir/phase_0.best.sorted.bam $dir/phase_0.best.bam
samtools sort -o $dir/phase_1.best.sorted.bam $dir/phase_1.best.bam
rm $dir/phase_*.sam
rm $dir/*best.bam
samtools index $dir/phase_0.best.sorted.bam
samtools index $dir/phase_1.best.sorted.bam

canu -pacbio-raw $dir/phase_0.fq -p A -d $dir/phase_0 genomeSize=8k useGrid=false
canu -pacbio-raw $dir/phase_1.fq -p A -d $dir/phase_1 genomeSize=8k useGrid=false
bwa mem $refe $dir/phase_0/A.contigs.fasta >  $dir/phase_0.contigs.sam
bwa mem $refe $dir/phase_1/A.contigs.fasta >  $dir/phase_1.contigs.sam
$trans $dir/phase_0.contigs.sam > $dir/phase_0.contigs.bam
$trans $dir/phase_1.contigs.sam > $dir/phase_1.contigs.bam
rm $dir/phase_0.contigs.sam
rm $dir/phase_1.contigs.sam
samtools index $dir/phase_0.contigs.bam
samtools index $dir/phase_1.contigs.bam

cp $dir/*.bam ~/D*
cp $dir/*.bai ~/D*
