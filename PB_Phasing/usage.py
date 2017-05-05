#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' Help information of the program'''

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def usage():
    ''' Help and usage information. '''
    info = '''
Program: phasing.py 
Version: 0.2.0
Contact: Zhou Ze <zhouze@genomocis.cn>

usage: phasing.py -b *mergedCCS.bestHitted.sorted.bam -o /output/directory -s start_coor -e end_coor
\t-h For help information.
\t-b Input BAM file, only sorted by position is allowed (samtools sort default paramter).
\t-o Output directory.
\t-m Target chromesome Name in BAM/SAM file.(default = \'chr6\')
\t-s Start coordinate of the region.
\t-e End coordinate of the region.
\t-p Upper heter-SNP proportion cutoff, alt fre/seq depth.(default = 0.75)
\t-d Lowwer heter-SNP proportion cutoff, alt fre/seq depth.(default = 0.25)
    '''
    print(info)

def check(input, output, chrom, reg_s, reg_e, max_heter, min_heter):
    ''' Paramters check.'''
    import sys
    if reg_e <= reg_s:
        print('\n*** Error -e value should > -s value')
        usage()
        sys.exit()
    elif min_heter > max_heter or max_heter > 1:
        print('\n*** Error: -p value should > -d value')
        usage()
        sys.exit()
    return 0

def author():
    ''' Return the author's name in Chinese.'''
    print('''Auther:\n\t_________  _________\n\t| __|__ |  |___|___|\n\t|___|___|  |___|___|\n\t\
|  ___  |      |      \n\t| |___| |      |      \n\t/       |      |     ''')
    
if __name__ == '__main__':
    check('in', 'out', 'chr6', 100, 1000, 0.5, 0.6)
