# Binary Tree based diploid genome Phasing software (BTP)

- PacBio CCS reads haplotype phasing program.
- Only for diploid genome (e.g. human genome).
- It performs better than SAMtools-phase on PacBio ultra long and high-error reads.

## Install:
```
git clone git@github.com:dazhouze/BTP.git
```

## Require:
Python3 and pysam module.

## Run:
- Firstly, generate accurate circular consensus sequence (CCS) from raw PacBio .h5 reads.
- Secondly, map CCS to reference genome through BWA-MEM (-x pacbio) or any other aligner. 
- Thirdly, select the highest AS score alignment result of read, if there are multiple-alignment results of same read.
- Fourthly, convert alignment results, SAM files, to BAM format, then merge, sort and bulid index for the BAM file.
- Fifthly, phase haplotypes by running BTP using command below, user have to assign a specific region of reference genome:
```
/path/to/phasing.py -m <chromosome_id> -s <start_pos> -e <end_pos> -o <output_dir> -b <aln.srt.bam>
```
For more information please see /path/to/phasing.py -h.

Example srcipt can be found in the example dirctory of this repo.

## Result:
Phased read Qname will be in your output_dir folder. Then you need to fetch FASTQ data by the Qname and assembly the FASTQ data by [Canu](https://github.com/marbl/canu) or any other assembler. The log files will be in your output_dir/log folder.

## TODO:
In multiple aligned reads, we select the alignment in highest AS score to avoid mis-alignment influencing phasing process. And in structure variants, there will be two true split alignments of a single reads, and in this case we will ignore one of the true alignment and generate a break point between phased blocks. Further, we will try to support this case.
