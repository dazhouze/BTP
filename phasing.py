#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.0.1'

'''
Python version phasing program.
More rebost.
'''


# Data structure 
'''
    PL = PositionalList()
    for x in ' WORLD':
        p = PL.add_after(p, x)
    print(PL.first().getNode().getElement())
    for x in PL:
        print(x, end = ' ')
    pode def init(self, qname, start, end, snp_pos=None, snp_alt=None, snp_qual=None, nextN):
'''
# Detect all SNP in all read 
# Identify heterozygous SNP marker, seq error and homo SNP 
# Detect ref-allel SNP in all read 
# Set seed 
# Initialize phase_0 SNP 
# Grow SNP tree (phase 0 SNP markers) 
# Filter phase_0 heter SNP markers 
# Scoring all reads 
# Determine 2 haplotype 
# Veen check 

class Read(object):
    #__slots__ = '__qname', '__start', '__end', '__snp_pos', '__snp_alt', '__snp_qual' #streamline memeory usage
    def __init__(self, qname, start, end, snp_pos=[None], snp_alt=[None], snp_qual=[None]):
        self.__qname = qname
        self.__start = start
        self.__end   = end
        self.__snp_pos  = snp_pos
        self.__snp_alt  = snp_alt
        self.__snp_qual = snp_qual

    def getQname(self):
        return self.__qname
    def getStart(self):
        return self.__start
    def getEnd(self):
        return self.__end
    def getSnpPos(self):
        return self.__snp_pos
    def getSnpAlt(self):
        return self.__snp_alt
    def getSnpQual(self):
        return self.__snp_qual

def main(input, output, chrom, reg_s, reg_e, seed_win, seed_cut,  max_heter, min_heter, snp_coin, score_cut):
    '''
    input bam file path
    reg_s region start coordinate
    reg_e region end coordinate
    '''
    import os, pysam
    from PB_Phasing import posList
    PL = posList.PositionalList() # initialize a positional list
    p = PL.add_first(Read('Begin', 0, 0)) # initialize the first item in PL list and store the position in varible p
    for x in PL:
        print(x.getQname())

    #Node def init(self, qname, start, end, snp_pos=None, snp_alt=None, snp_qual=None, nextN):

    # Data structure 
    # Detect all SNP in all read 
    ##### Identify SAM/BAM file to open different pysam IO handle. #####
    if os.path.splitext(input)[1] == '.sam': # SAM file
        bamfile = pysam.AlignmentFile(input, "r")
    elif os.path.splitext(input)[1] == '.bam': # BAM file
        bamfile = pysam.AlignmentFile(input, "rb") # file handle of BAM file
    else:
        raise IOError('Please choose a SAM/BAM file!')
    target = bamfile.fetch('chr6', reg_s, reg_e) # iterable method of target region read
    for read in target:
        for p in range(read.reference_start, read.reference_end):
            ind = p - reg_s # index
            if ind < reg_l:
               depAy[ind] += 1



def usage():
    info = '''
Program: phasing.py 
Version: 0.2.0
Contact: Zhou Ze <zhouze@genomocis.cn>

Usage: phasing.py -b *mergedCCS.bestHitted.sorted.bam -o /output/directory -s start_coor -e end_coor
\t-h For more information.
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
    import getopt, sys
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:o:c:u:w:p:d:v:s:e:m:")
    except getopt.GetoptError as err:
        usage()
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)

    output = None # output dirctory
    input = '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/merged5YH.best.ccs.sort.bam' # input BAM/SAM file
    temp_delete = False # if delete the temp direcory
    chrom = 'chr6' # chromesome name
    reg_s = 0 # start coordinate of the region
    reg_e = 0 # end coordinate of the region
    max_heter = 0.75 # upper heter snp cutoff, alt fre/seq depth
    min_heter = 0.24 # #lowwer heter snp cutoff, alt fre/seq depth
    seed_win = 300 # window size of seed region selection
    seed_cut = 0.30 # cutoff value of seed (SNP) pattern selection
    snp_coin = 0.90  # SNP coincide proportion when extending
    score_cut = 0.55 # phase 0 and 1 cutoff value of scoring

    for o, a in opts:
        if o == '-h':
            usage()
            sys.exit()
        elif o == '-b':
            input = a
        elif o == '-o':
            output = a

        elif o == '-m':
            chrom = a # chromesome
        elif o == '-s':
            reg_s = a # start coordinate of the region
        elif o == '-e':
            reg_e = a # end coordinate of the region

        elif o == '-p':
            max_heter = a # upper heter snp cutoff, alt fre/seq depth
        elif o == '-d':
            min_heter = a # #lowwer heter snp cutoff, alt fre/seq depth

        elif o == '-w':
            seed_win = a # window size of seed region selection
        elif o == '-v':
            seed_cut = a # cutoff value of seed (SNP) pattern selection

        elif o == '-c':
            snp_coin = a # SNP coincide proportion when extending
        elif o == '-u':
            score_cut = a # phase 0 and 1 cutoff value of scoring
    main(input, output, chrom, reg_s, reg_e, seed_win, seed_cut,  max_heter, min_heter, snp_coin, score_cut)
