#!/bin/bash
# low long region 91k
#python phasing.py -s 28497412   -e 28588898
dir=./test
rm -f $dir/*
# A gene
python phasing.py -s 29910247   -e 29913661 -o $dir > $dir/hit.txt
# 7 k
#python phasing.py -s 29053005 -e 29060759 -o $dir
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
samtools index $dir/phase_0.best.sorted.bam
samtools index $dir/phase_1.best.sorted.bam

cp $dir/phase_0.best.sorted.bam ~/D*
cp $dir/phase_1.best.sorted.bam ~/D*
cp $dir/phase_0.best.sorted.bam.bai ~/D*
cp $dir/phase_1.best.sorted.bam.bai ~/D*
cp $dir/hit.txt ~/D*
