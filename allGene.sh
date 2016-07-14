rm -r A
rm -r B
rm -r C
rm -r DQA1
rm -r DQB1
rm -r DPAB1
rm -r DRB1
perl perl.Phasing.pl  -o A ../../6_filtered_result/output/BAM/A.sort.bam
perl perl.qname2fastq.pl A/phase.0.qname
perl perl.qname2fastq.pl A/phase.1.qname
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa A/phase.0.qname.fastq > A/A.0.fastq.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa A/phase.1.qname.fastq > A/A.1.fastq.sam
samtools view -Sb A/A.1.fastq.sam | samtools sort > A/A.1.fastq.bam
samtools view -Sb A/A.0.fastq.sam | samtools sort > A/A.0.fastq.bam
samtools index A/A.0.fastq.bam
samtools index A/A.1.fastq.bam
canu  useGrid=false -p A -d A/A_0 genomesize=6000 -pacbio-raw  A/phase.0.qname.fastq
canu  useGrid=false -p A -d A/A_1 genomesize=6000 -pacbio-raw  A/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa A/A_0/A.consensus.fasta > A/A.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa A/A_1/A.consensus.fasta > A/A.1.contig.sam
samtools view -Sb A/A.1.contig.sam | samtools sort > A/A.1.contig.bam
samtools view -Sb A/A.0.contig.sam | samtools sort > A/A.0.contig.bam
samtools index A/A.0.contig.bam
samtools index A/A.1.contig.bam

perl perl.Phasing.pl  -o B ../../6_filtered_result/output/BAM/B.sort.bam
perl perl.qname2fastq.pl B/phase.0.qname
perl perl.qname2fastq.pl B/phase.1.qname
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa B/phase.0.qname.fastq > B/B.0.fastq.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa B/phase.1.qname.fastq > B/B.1.fastq.sam
samtools view -Sb B/B.1.fastq.sam | samtools sort > B/B.1.fastq.bam
samtools view -Sb B/B.0.fastq.sam | samtools sort > B/B.0.fastq.bam
samtools index B/B.0.fastq.bam
samtools index B/B.1.fastq.bam
canu  useGrid=false -p B -d B/B_0 genomesize=6000 -pacbio-raw  B/phase.0.qname.fastq
canu  useGrid=false -p B -d B/B_1 genomesize=6000 -pacbio-raw  B/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa B/B_0/B.consensus.fasta > B/B.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa B/B_1/B.consensus.fasta > B/B.1.contig.sam
samtools view -Sb B/B.1.contig.sam | samtools sort > B/B.1.contig.bam
samtools view -Sb B/B.0.contig.sam | samtools sort > B/B.0.contig.bam
samtools index B/B.0.contig.bam
samtools index B/B.1.contig.bam

perl perl.Phasing.pl  -o C ../../6_filtered_result/output/BAM/C.sort.bam
perl perl.qname2fastq.pl C/phase.0.qname
perl perl.qname2fastq.pl C/phase.1.qname
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa C/phase.0.qname.fastq > C/C.0.fastq.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa C/phase.1.qname.fastq > C/C.1.fastq.sam
samtools view -Sb C/C.1.fastq.sam | samtools sort > C/C.1.fastq.bam
samtools view -Sb C/C.0.fastq.sam | samtools sort > C/C.0.fastq.bam
samtools index C/C.0.fastq.bam
samtools index C/C.1.fastq.bam
canu  useGrid=false -p C -d C/C_0 genomesize=6000 -pacbio-raw  C/phase.0.qname.fastq
canu  useGrid=false -p C -d C/C_1 genomesize=6000 -pacbio-raw  C/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa C/C_0/C.consensus.fasta > C/C.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa C/C_1/C.consensus.fasta > C/C.1.contig.sam
samtools view -Sb C/C.1.contig.sam | samtools sort > C/C.1.contig.bam
samtools view -Sb C/C.0.contig.sam | samtools sort > C/C.0.contig.bam
samtools index C/C.0.contig.bam
samtools index C/C.1.contig.bam

perl perl.Phasing.pl  -o DQA1 ../../6_filtered_result/output/BAM/DQA1.sort.bam
perl perl.qname2fastq.pl DQA1/phase.0.qname
perl perl.qname2fastq.pl DQA1/phase.1.qname
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQA1/phase.0.qname.fastq > DQA1/DQA1.0.fastq.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQA1/phase.1.qname.fastq > DQA1/DQA1.1.fastq.sam
samtools view -Sb DQA1/DQA1.1.fastq.sam | samtools sort > DQA1/DQA1.1.fastq.bam
samtools view -Sb DQA1/DQA1.0.fastq.sam | samtools sort > DQA1/DQA1.0.fastq.bam
samtools index DQA1/DQA1.0.fastq.bam
samtools index DQA1/DQA1.1.fastq.bam
canu  useGrid=false -p DQA1 -d DQA1/DQA1_0 genomesize=6000 -pacbio-raw  DQA1/phase.0.qname.fastq
canu  useGrid=false -p DQA1 -d DQA1/DQA1_1 genomesize=6000 -pacbio-raw  DQA1/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQA1/DQA1_0/DQA1.consensus.fasta > DQA1/DQA1.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQA1/DQA1_1/DQA1.consensus.fasta > DQA1/DQA1.1.contig.sam
samtools view -Sb DQA1/DQA1.1.contig.sam | samtools sort > DQA1/DQA1.1.contig.bam
samtools view -Sb DQA1/DQA1.0.contig.sam | samtools sort > DQA1/DQA1.0.contig.bam
samtools index DQA1/DQA1.0.contig.bam
samtools index DQA1/DQA1.1.contig.bam

perl perl.Phasing.pl  -o DQB1 ../../6_filtered_result/output/BAM/DQB1.sort.bam
perl perl.qname2fastq.pl DQB1/phase.0.qname
perl perl.qname2fastq.pl DQB1/phase.1.qname
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQB1/phase.0.qname.fastq > DQB1/DQB1.0.fastq.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQB1/phase.1.qname.fastq > DQB1/DQB1.1.fastq.sam
samtools view -Sb DQB1/DQB1.1.fastq.sam | samtools sort > DQB1/DQB1.1.fastq.bam
samtools view -Sb DQB1/DQB1.0.fastq.sam | samtools sort > DQB1/DQB1.0.fastq.bam
samtools index DQB1/DQB1.0.fastq.bam
samtools index DQB1/DQB1.1.fastq.bam
canu  useGrid=false -p DQB1 -d DQB1/DQB1_0 genomesize=6000 -pacbio-raw  DQB1/phase.0.qname.fastq
canu  useGrid=false -p DQB1 -d DQB1/DQB1_1 genomesize=6000 -pacbio-raw  DQB1/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQB1/DQB1_0/DQB1.consensus.fasta > DQB1/DQB1.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DQB1/DQB1_1/DQB1.consensus.fasta > DQB1/DQB1.1.contig.sam
samtools view -Sb DQB1/DQB1.1.contig.sam | samtools sort > DQB1/DQB1.1.contig.bam
samtools view -Sb DQB1/DQB1.0.contig.sam | samtools sort > DQB1/DQB1.0.contig.bam
samtools index DQB1/DQB1.0.contig.bam
samtools index DQB1/DQB1.1.contig.bam

perl perl.Phasing.pl  -o DPAB1 ../../6_filtered_result/output/BAM/DPAB1.sort.bam
perl perl.qname2fastq.pl DPAB1/phase.0.qname
perl perl.qname2fastq.pl DPAB1/phase.1.qname
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DPAB1/phase.0.qname.fastq > DPAB1/DPAB1.0.fastq.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DPAB1/phase.1.qname.fastq > DPAB1/DPAB1.1.fastq.sam
samtools view -Sb DPAB1/DPAB1.1.fastq.sam | samtools sort > DPAB1/DPAB1.1.fastq.bam
samtools view -Sb DPAB1/DPAB1.0.fastq.sam | samtools sort > DPAB1/DPAB1.0.fastq.bam
samtools index DPAB1/DPAB1.0.fastq.bam
samtools index DPAB1/DPAB1.1.fastq.bam
canu  useGrid=false -p DPAB1 -d DPAB1/DPAB1_0 genomesize=6000 -pacbio-raw  DPAB1/phase.0.qname.fastq
canu  useGrid=false -p DPAB1 -d DPAB1/DPAB1_1 genomesize=6000 -pacbio-raw  DPAB1/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DPAB1/DPAB1_0/DPAB1.consensus.fasta > DPAB1/DPAB1.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DPAB1/DPAB1_1/DPAB1.consensus.fasta > DPAB1/DPAB1.1.contig.sam
samtools view -Sb DPAB1/DPAB1.1.contig.sam | samtools sort > DPAB1/DPAB1.1.contig.bam
samtools view -Sb DPAB1/DPAB1.0.contig.sam | samtools sort > DPAB1/DPAB1.0.contig.bam
samtools index DPAB1/DPAB1.0.contig.bam
samtools index DPAB1/DPAB1.1.contig.bam

perl perl.Phasing.pl  -o DRB1 ../../6_filtered_result/output/BAM/DRB1.sort.bam
perl perl.qname2fastq.pl DRB1/phase.0.qname
perl perl.qname2fastq.pl DRB1/phase.1.qname
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DRB1/phase.0.qname.fastq > DRB1/DRB1.0.fastq.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DRB1/phase.1.qname.fastq > DRB1/DRB1.1.fastq.sam
samtools view -Sb DRB1/DRB1.1.fastq.sam | samtools sort > DRB1/DRB1.1.fastq.bam
samtools view -Sb DRB1/DRB1.0.fastq.sam | samtools sort > DRB1/DRB1.0.fastq.bam
samtools index DRB1/DRB1.0.fastq.bam
samtools index DRB1/DRB1.1.fastq.bam
canu  useGrid=false -p DRB1 -d DRB1/DRB1_0 genomesize=6000 -pacbio-raw  DRB1/phase.0.qname.fastq
canu  useGrid=false -p DRB1 -d DRB1/DRB1_1 genomesize=6000 -pacbio-raw  DRB1/phase.1.qname.fastq
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DRB1/DRB1_0/DRB1.consensus.fasta > DRB1/DRB1.0.contig.sam
bwa mem ~/team/dataBase/hg19_chr6/chr6.fa DRB1/DRB1_1/DRB1.consensus.fasta > DRB1/DRB1.1.contig.sam
samtools view -Sb DRB1/DRB1.1.contig.sam | samtools sort > DRB1/DRB1.1.contig.bam
samtools view -Sb DRB1/DRB1.0.contig.sam | samtools sort > DRB1/DRB1.0.contig.bam
samtools index DRB1/DRB1.0.contig.bam
samtools index DRB1/DRB1.1.contig.bam
