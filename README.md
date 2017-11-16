# Binary tree based diploid genome Phasing software

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
- Firstly, generate accurate circular consensus sequenct (CCS) from raw PacBio .h5 reads.
- Secondly, map CCS to reference genome through BWA-MEM (-x pacbio) or any other aligner.
- Thirdly, convert alignment files to BAM files, then merge, sort and bulid index of BAM file.
- Fourthly, run BTP using command below:
```
/path/to/phasing.py -m chromosome -s start_pos -e end_pos -o output_dir -b BAM
```
For more information please see /path/to/phasing.py -h
Example srcipt can be found in the example dirctory of this repo.

## Result:
Phased read Qname will be in your output_dir folder. Then you need to fetch FASTQ data by the Qname and assembly the FASTQ data by [Canu](https://github.com/marbl/canu) or any other assembler. The log files will be in your output_dir/log folder.
