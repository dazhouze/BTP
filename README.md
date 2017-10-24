# Binary tree based diploid genome Phasing

- PacBio CCS reads haplotype phasing program.
- Perform better than Samtools phase on PacBio extra long and high error rate reads.

## Install
```
git clone git@github.com:dazhouze/BTP.git
```

## Require
Python3 and pysam module

## Run
Firstly, get ccs corrected Fastq files from raw PacBio reads.
Secondly, get sam files by running BWA-MEM (-x pacbio) alignment to ccs Fastq files.
Thirdly, convert sam files to bam files, merge them, sort merged bam file and index it.
Then run this program by:
```
/path/to/phasing.py -m chromosome -s start_pos -e end_pos -o output_dir
```
For more information: /path/to/phasing.py -h
Example srcipt can be found in the example dirctory of this repo.

## Result
Phased fragments' QNAME will be in your output_dir folder. Then you need to fetch FASTQ file 
Logs file will be in your output_dir/log folder.
