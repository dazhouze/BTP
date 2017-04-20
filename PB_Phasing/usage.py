#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def Usage():
    info = '''
Program: phasing.py 
Version: 0.2.0
Contact: Zhou Ze <zhouze@genomocis.cn>

Usage: phasing.py -b *mergedCCS.bestHitted.sorted.bam -o /output/directory -s start_coor -e end_coor
\t-h For help information.
\t-b Input BAM file, only sorted by position is allowed (samtools sort default paramter).
\t-o Output directory.
\t-m Target chromesome Name in BAM/SAM file.(default = \'chr6\')
\t-s Start coordinate of the region.
\t-e End coordinate of the region.
\t-p Upper heter-SNP proportion cutoff, alt fre/seq depth.(default = 0.75)
\t-d Lowwer heter-SNP proportion cutoff, alt fre/seq depth.(default = 0.25)
\t-u Phase 0 and 1 > score cutoff. 0.0--1.0(score_cut = 0.55)
    '''
    print(info)
def check(input, output, chrom, reg_s, reg_e, max_heter, min_heter, score_cut):
    import sys
    if reg_e - reg_s <= 2000:
        Usage()
        sys.exit()
    elif min_heter > max_heter > 1:
        Usage()
        sys.exit()
    return 0

if __name__ == '__main__':
    Usage()
