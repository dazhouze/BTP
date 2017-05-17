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
\t-r Sequencing error SNP proportion cutoff, alt_allele/seq_depth.(default = 0.25)
    '''
    print(info)

def check(input, output, chrom, reg_s, reg_e, seq_error):
    ''' Paramters check.'''
    import sys
    if reg_e <= reg_s:
        print('\n*** Parameter error: -e value should > -s value')
        usage()
        sys.exit()
    elif seq_error < 0:
        print('\n*** Parameter error: -r value should > 0')
        usage()
        sys.exit()
    elif seq_error > 0.5:
        print('\n*** Parameter error: -r value should < 0.5')
        usage()
        sys.exit()
    return 0

def author():
    ''' Return the author's name in Chinese.'''
    print('''Auther:\n\t_________  _________\n\t| __|__ |  |___|___|\n\t|___|___|  |___|___|\n\t\
|  ___  |      |      \n\t| |___| |      |      \n\t/       |      |     ''')
    
if __name__ == '__main__':
    check('in', 'out', 'chr6', 100, 1000, 0.3)
