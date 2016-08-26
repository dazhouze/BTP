rm -rf DRB1
#./Phasing -o DRB1 /home/zhouze/team/YH05_PacBio/YH06_YH07/2_bestHit/YH06.mem-x.DRB1.bam
perl ./perl.Phasing.pl -o DRB1 -w 400 -e 0.4 /home/zhouze/team/YH05_PacBio/YH06_YH07/2_bestHit/YH06.mem-x.DRB1.bam
python Qname2Fastq.py DRB1/phase.0.qname ~/YH06_YH07/1_SMRT/YH06/data/reads_of_insert.fastq > DRB1/phase.0.qname.fastq
python Qname2Fastq.py DRB1/phase.1.qname ~/YH06_YH07/1_SMRT/YH06/data/reads_of_insert.fastq > DRB1/phase.1.qname.fastq
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
