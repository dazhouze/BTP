#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'Zhou Ze'
__version__ = '0.0.1'

def usage():
    info = '''
Program: phasing.py 
Version: 0.2.0
Contact: Zhou Ze <zhouze@genomocis.cn>

Usage: phasing.py -b *mergedCCS.bestHitted.sorted.bam -o /output/directory -s start_coor -e end_coor
\t-h For help information.
\t-b Input BAM/SAM file.
\t-o Output directory.
\t-t Delete the TEMP directory.
\t-m Target chromesome Name in BAM/SAM file.(chrom = \'chr6\')
\t-s Start coordinate of the region.
\t-e End coordinate of the region.
\t-p Upper heter-SNP proportion cutoff, alt fre/seq depth.(max_heter = 0.75)
\t-d Lowwer heter-SNP proportion cutoff, alt fre/seq depth.(min_heter = 0.24)
\t-w Window size of seed selection.(seed_win = 300)
\t-v Seed (SNP) pattern selection cutoff value.(seed_cut = 0.30)
\t-c Coincident SNP proportion. 0.0-1.0(snp_coin = 0.90)
\t-u Phase 0 and 1 > score cutoff. 0.0--1.0(score_cut = 0.55)
    '''
    print(info)

if __name__ == '__main__':
    usage()
