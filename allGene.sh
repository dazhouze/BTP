#!/bin/bash

#rm previous directory
rm -rf A
rm -rf B
rm -rf C
rm -rf DQA1
rm -rf DQB1
rm -rf DPA1
rm -rf DPB1
rm -rf DRB1


#Gene region phase
perl perl.Phasing.pl  -o A -c 0.92 ../../6_filtered_result/output/BAM/A.sort.bam
perl perl.Phasing.pl  -o B ../../6_filtered_result/output/BAM/B.sort.bam
perl perl.Phasing.pl  -o C ../../6_filtered_result/output/BAM/C.sort.bam
perl perl.Phasing.pl  -o DQA1 -w 300 ../../6_filtered_result/output/BAM/DQA1.sort.bam
perl perl.Phasing.pl  -o DQB1 -w 300 ../../6_filtered_result/output/BAM/DQB1.sort.bam
perl perl.Phasing.pl  -o DPA1 -c 0.86 -w 300 ../../6_filtered_result/output/BAM/DPA1.sort.bam
perl perl.Phasing.pl  -o DPB1 -w 300 ../../6_filtered_result/output/BAM/DPB1.sort.bam
perl perl.Phasing.pl  -o DRB1 ../../6_filtered_result/output/BAM/DRB1.sort.bam

exit 0

#Gene region qname to fastq and alignment
##Gene A
perl perl.qname2fastq.pl A/phase.0.qname
perl perl.qname2fastq.pl A/phase.1.qname
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa A/phase.0.qname.fastq > A/A.0.fastq.sam
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa A/phase.1.qname.fastq > A/A.1.fastq.sam
samtools view -Sb A/A.1.fastq.sam | samtools sort > A/A.1.fastq.bam
samtools view -Sb A/A.0.fastq.sam | samtools sort > A/A.0.fastq.bam
samtools index A/A.0.fastq.bam
samtools index A/A.1.fastq.bam
##Gene B
perl perl.qname2fastq.pl B/phase.0.qname
perl perl.qname2fastq.pl B/phase.1.qname
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa B/phase.0.qname.fastq > B/B.0.fastq.sam
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa B/phase.1.qname.fastq > B/B.1.fastq.sam
samtools view -Sb B/B.1.fastq.sam | samtools sort > B/B.1.fastq.bam
samtools view -Sb B/B.0.fastq.sam | samtools sort > B/B.0.fastq.bam
samtools index B/B.0.fastq.bam
samtools index B/B.1.fastq.bam
##Gene C
perl perl.qname2fastq.pl C/phase.0.qname
perl perl.qname2fastq.pl C/phase.1.qname
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa C/phase.0.qname.fastq > C/C.0.fastq.sam
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa C/phase.1.qname.fastq > C/C.1.fastq.sam
samtools view -Sb C/C.1.fastq.sam | samtools sort > C/C.1.fastq.bam
samtools view -Sb C/C.0.fastq.sam | samtools sort > C/C.0.fastq.bam
samtools index C/C.0.fastq.bam
samtools index C/C.1.fastq.bam
##Gene DQA1
perl perl.qname2fastq.pl DQA1/phase.0.qname
perl perl.qname2fastq.pl DQA1/phase.1.qname
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DQA1/phase.0.qname.fastq > DQA1/DQA1.0.fastq.sam
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DQA1/phase.1.qname.fastq > DQA1/DQA1.1.fastq.sam
samtools view -Sb DQA1/DQA1.1.fastq.sam | samtools sort > DQA1/DQA1.1.fastq.bam
samtools view -Sb DQA1/DQA1.0.fastq.sam | samtools sort > DQA1/DQA1.0.fastq.bam
samtools index DQA1/DQA1.0.fastq.bam
samtools index DQA1/DQA1.1.fastq.bam
##Gene DQB1
perl perl.qname2fastq.pl DQB1/phase.0.qname
perl perl.qname2fastq.pl DQB1/phase.1.qname
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DQB1/phase.0.qname.fastq > DQB1/DQB1.0.fastq.sam
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DQB1/phase.1.qname.fastq > DQB1/DQB1.1.fastq.sam
samtools view -Sb DQB1/DQB1.1.fastq.sam | samtools sort > DQB1/DQB1.1.fastq.bam
samtools view -Sb DQB1/DQB1.0.fastq.sam | samtools sort > DQB1/DQB1.0.fastq.bam
samtools index DQB1/DQB1.0.fastq.bam
samtools index DQB1/DQB1.1.fastq.bam
##Gene DPA1
perl perl.qname2fastq.pl DPA1/phase.0.qname
perl perl.qname2fastq.pl DPA1/phase.1.qname
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DPA1/phase.0.qname.fastq > DPA1/DPA1.0.fastq.sam
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DPA1/phase.1.qname.fastq > DPA1/DPA1.1.fastq.sam
samtools view -Sb DPA1/DPA1.1.fastq.sam | samtools sort > DPA1/DPA1.1.fastq.bam
samtools view -Sb DPA1/DPA1.0.fastq.sam | samtools sort > DPA1/DPA1.0.fastq.bam
samtools index DPA1/DPA1.0.fastq.bam
samtools index DPA1/DPA1.1.fastq.bam
##Gene DPB1
perl perl.qname2fastq.pl DPB1/phase.0.qname
perl perl.qname2fastq.pl DPB1/phase.1.qname
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DPB1/phase.0.qname.fastq > DPB1/DPB1.0.fastq.sam
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DPB1/phase.1.qname.fastq > DPB1/DPB1.1.fastq.sam
samtools view -Sb DPB1/DPB1.1.fastq.sam | samtools sort > DPB1/DPB1.1.fastq.bam
samtools view -Sb DPB1/DPB1.0.fastq.sam | samtools sort > DPB1/DPB1.0.fastq.bam
samtools index DPB1/DPB1.0.fastq.bam
samtools index DPB1/DPB1.1.fastq.bam
##Gene DRB1
perl perl.qname2fastq.pl DRB1/phase.0.qname
perl perl.qname2fastq.pl DRB1/phase.1.qname
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DRB1/phase.0.qname.fastq > DRB1/DRB1.0.fastq.sam
bwa mem -x pacbio ~/team/dataBase/hg19_chr6/chr6.fa DRB1/phase.1.qname.fastq > DRB1/DRB1.1.fastq.sam
samtools view -Sb DRB1/DRB1.1.fastq.sam | samtools sort > DRB1/DRB1.1.fastq.bam
samtools view -Sb DRB1/DRB1.0.fastq.sam | samtools sort > DRB1/DRB1.0.fastq.bam
samtools index DRB1/DRB1.0.fastq.bam
samtools index DRB1/DRB1.1.fastq.bam

#Two phase Canu assembly
##Gene A
canu  useGrid=false -p A -d A/A_0 genomesize=6000 -pacbio-raw  A/phase.0.qname.fastq
canu  useGrid=false -p A -d A/A_1 genomesize=6000 -pacbio-raw  A/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa A/A_0/A.consensus.fasta > A/A.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa A/A_1/A.consensus.fasta > A/A.1.contig.sam
samtools view -Sb A/A.1.contig.sam | samtools sort > A/A.1.contig.bam
samtools view -Sb A/A.0.contig.sam | samtools sort > A/A.0.contig.bam
samtools index A/A.0.contig.bam
samtools index A/A.1.contig.bam
##Gene B
canu  useGrid=false -p B -d B/B_0 genomesize=6000 -pacbio-raw  B/phase.0.qname.fastq
canu  useGrid=false -p B -d B/B_1 genomesize=6000 -pacbio-raw  B/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa B/B_0/B.consensus.fasta > B/B.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa B/B_1/B.consensus.fasta > B/B.1.contig.sam
samtools view -Sb B/B.1.contig.sam | samtools sort > B/B.1.contig.bam
samtools view -Sb B/B.0.contig.sam | samtools sort > B/B.0.contig.bam
samtools index B/B.0.contig.bam
samtools index B/B.1.contig.bam
##Gene C
canu  useGrid=false -p C -d C/C_0 genomesize=6000 -pacbio-raw  C/phase.0.qname.fastq
canu  useGrid=false -p C -d C/C_1 genomesize=6000 -pacbio-raw  C/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa C/C_0/C.consensus.fasta > C/C.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa C/C_1/C.consensus.fasta > C/C.1.contig.sam
samtools view -Sb C/C.1.contig.sam | samtools sort > C/C.1.contig.bam
samtools view -Sb C/C.0.contig.sam | samtools sort > C/C.0.contig.bam
samtools index C/C.0.contig.bam
samtools index C/C.1.contig.bam
##Gene DQA1
canu  useGrid=false -p DQA1 -d DQA1/DQA1_0 genomesize=6000 -pacbio-raw  DQA1/phase.0.qname.fastq
canu  useGrid=false -p DQA1 -d DQA1/DQA1_1 genomesize=6000 -pacbio-raw  DQA1/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQA1/DQA1_0/DQA1.consensus.fasta > DQA1/DQA1.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQA1/DQA1_1/DQA1.consensus.fasta > DQA1/DQA1.1.contig.sam
samtools view -Sb DQA1/DQA1.1.contig.sam | samtools sort > DQA1/DQA1.1.contig.bam
samtools view -Sb DQA1/DQA1.0.contig.sam | samtools sort > DQA1/DQA1.0.contig.bam
samtools index DQA1/DQA1.0.contig.bam
samtools index DQA1/DQA1.1.contig.bam
##Gene DQB1
canu  useGrid=false -p DQB1 -d DQB1/DQB1_0 genomesize=6000 -pacbio-raw  DQB1/phase.0.qname.fastq
canu  useGrid=false -p DQB1 -d DQB1/DQB1_1 genomesize=6000 -pacbio-raw  DQB1/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQB1/DQB1_0/DQB1.consensus.fasta > DQB1/DQB1.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQB1/DQB1_1/DQB1.consensus.fasta > DQB1/DQB1.1.contig.sam
samtools view -Sb DQB1/DQB1.1.contig.sam | samtools sort > DQB1/DQB1.1.contig.bam
samtools view -Sb DQB1/DQB1.0.contig.sam | samtools sort > DQB1/DQB1.0.contig.bam
samtools index DQB1/DQB1.0.contig.bam
samtools index DQB1/DQB1.1.contig.bam
##Gene DPA1
canu  useGrid=false -p DPA1 -d DPA1/DPA1_0 genomesize=6000 -pacbio-raw  DPA1/phase.0.qname.fastq
canu  useGrid=false -p DPA1 -d DPA1/DPA1_1 genomesize=6000 -pacbio-raw  DPA1/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DPA1/DPA1_0/DPA1.consensus.fasta > DPA1/DPA1.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DPA1/DPA1_1/DPA1.consensus.fasta > DPA1/DPA1.1.contig.sam
samtools view -Sb DPA1/DPA1.1.contig.sam | samtools sort > DPA1/DPA1.1.contig.bam
samtools view -Sb DPA1/DPA1.0.contig.sam | samtools sort > DPA1/DPA1.0.contig.bam
samtools index DPA1/DPA1.0.contig.bam
samtools index DPA1/DPA1.1.contig.bam
##Gene DPB1
canu  useGrid=false -p DPB1 -d DPB1/DPB1_0 genomesize=6000 -pacbio-raw  DPB1/phase.0.qname.fastq
canu  useGrid=false -p DPB1 -d DPB1/DPB1_1 genomesize=6000 -pacbio-raw  DPB1/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DPB1/DPB1_0/DPB1.consensus.fasta > DPB1/DPB1.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DPB1/DPB1_1/DPB1.consensus.fasta > DPB1/DPB1.1.contig.sam
samtools view -Sb DPB1/DPB1.1.contig.sam | samtools sort > DPB1/DPB1.1.contig.bam
samtools view -Sb DPB1/DPB1.0.contig.sam | samtools sort > DPB1/DPB1.0.contig.bam
samtools index DPB1/DPB1.0.contig.bam
samtools index DPB1/DPB1.1.contig.bam
##Gene DRB1
canu  useGrid=false -p DRB1 -d DRB1/DRB1_0 genomesize=6000 -pacbio-raw  DRB1/phase.0.qname.fastq
canu  useGrid=false -p DRB1 -d DRB1/DRB1_1 genomesize=6000 -pacbio-raw  DRB1/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DRB1/DRB1_0/DRB1.consensus.fasta > DRB1/DRB1.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DRB1/DRB1_1/DRB1.consensus.fasta > DRB1/DRB1.1.contig.sam
samtools view -Sb DRB1/DRB1.1.contig.sam | samtools sort > DRB1/DRB1.1.contig.bam
samtools view -Sb DRB1/DRB1.0.contig.sam | samtools sort > DRB1/DRB1.0.contig.bam
samtools index DRB1/DRB1.0.contig.bam
samtools index DRB1/DRB1.1.contig.bam
