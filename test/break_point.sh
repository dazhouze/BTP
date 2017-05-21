dir=./Break_Point
gene=bp
rm -rf $dir
python3 phasing.py -s 29684448 -e 29689116 -o $dir 
wc -l $dir/*txt
# fetch fastq
python3 fetch.py $dir/phase_0.0.txt $dir/phase_1.0.txt $dir/phase_0.0.fq $dir/phase_1.0.fq
python3 fetch.py $dir/phase_0.1.txt $dir/phase_1.1.txt $dir/phase_0.1.fq $dir/phase_1.1.fq
wc -l $dir/*fq
# BWA Alignment Command
align="bwa mem -x pacbio"
# Homo Sapiens Reference Renome hg19
refe=/ifs1/ST_IM/USER/database/Resequence/database/Ref/hg19/ucsc.hg19.fa
# SAM File Best Hit Perl Script
best="/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/workflow/step_2_CCS2BAM/SAMbestHit.pl" 
# SAM To BAM Format Conversion Command
trans="samtools view -Sb"

$align $refe $dir/phase_0.0.fq > $dir/phase_0.0.sam
$align $refe $dir/phase_1.0.fq > $dir/phase_1.0.sam
perl $best $dir/phase_0.0.sam $dir/phase_0.0.best.sam
perl $best $dir/phase_1.0.sam $dir/phase_1.0.best.sam
$trans $dir/phase_0.0.best.sam > $dir/phase_0.0.best.bam
$trans $dir/phase_1.0.best.sam > $dir/phase_1.0.best.bam
samtools sort -o $dir/phase_0.0.best.sorted.bam $dir/phase_0.0.best.bam
samtools sort -o $dir/phase_1.0.best.sorted.bam $dir/phase_1.0.best.bam
rm $dir/phase_*.sam
rm $dir/*best.bam
samtools index $dir/phase_0.0.best.sorted.bam
samtools index $dir/phase_1.0.best.sorted.bam

$align $refe $dir/phase_0.1.fq > $dir/phase_0.1.sam
$align $refe $dir/phase_1.1.fq > $dir/phase_1.1.sam
perl $best $dir/phase_0.1.sam $dir/phase_0.1.best.sam
perl $best $dir/phase_1.1.sam $dir/phase_1.1.best.sam
$trans $dir/phase_0.1.best.sam > $dir/phase_0.1.best.bam
$trans $dir/phase_1.1.best.sam > $dir/phase_1.1.best.bam
samtools sort -o $dir/phase_0.1.best.sorted.bam $dir/phase_0.1.best.bam
samtools sort -o $dir/phase_1.1.best.sorted.bam $dir/phase_1.1.best.bam
rm $dir/phase_*.sam
rm $dir/*best.bam
samtools index $dir/phase_0.1.best.sorted.bam
samtools index $dir/phase_1.1.best.sorted.bam
