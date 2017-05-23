# PB_Phasing

- PacBio reads haplotype phasing program by applying binary tree.
- Perform better than Samtools phase on PacBio extra long reads.

## Install
```
git clone
```

## Run
Firstly, get ccs corrected Fastq files from raw PacBio reads.
Secondly, get sam files by running BWA-MEM (-x pacbio) alignment to ccs Fastq files.
Thirdly, convert sam files to bam files, merge them, sort merged bam file and index it.
Then run this program by:
```
/path/to/phasing.py -m chromosome -s start_pos -e end_pos -o output_dir
```
See more information by: -h

## Result
Phased fragments' qname will be in your output_dir folder.
Logs file will be in your output_dir/log folder.