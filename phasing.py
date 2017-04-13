#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'Zhou Ze'
__version__ = '0.0.1'
'''
Python version phasing program.
More rebost.
'''
# 1.Detect all SNP in all read 
# 2.Identify heterozygous SNP marker, seq error and homo SNP 
# 3.Detect ref-allel SNP in all read 
# 4.Set seed 
# 5.Initialize phase_0 SNP 
# 6.Grow SNP tree (phase 0 SNP markers) 
# 7.Filter phase_0 heter SNP markers 
# 8.Scoring all reads 
# 9.Determine 2 haplotype 
# 10.Veen check 

##### Data structure #####
class Read(object):
    __slots__ = '__qname', '__start', '__end', '__ave_sq', '__snp' #streamline memeory usage
    def __init__(self, qname, start, end, ave_sq, snp):
        self.__qname = qname # QNAME of read
        self.__start = start # start coordinate in reference genome
        self.__end   = end # end coordinate in reference genome
        self.__ave_sq  = ave_sq # average value of sequencing qulity
        self.__snp  = snp # dict of snp. key is position, vaule is [snp_alt, snp_qual]
    def getQname(self):
        return self.__qname
    def getStart(self):
        return self.__start
    def getEnd(self):
        return self.__end
    def getAveSq(self):
        return self.__ave_sq
    def getSnp(self):
        return self.__snp

def main(input, output, chrom, reg_s, reg_e, seed_win, seed_cut,  max_heter, min_heter, snp_coin, score_cut):
    '''
    input: SAM/BAM file path
    output: output directory path
    reg_s: region start coordinate
    reg_e: region end coordinate
    '''
    import os, pysam
    from PB_Phasing import posList, usasge, heterSnp

    # init varibls for Phasing #
    heter_snp = {} 
    ''' All possible position of SNPs.
    key is pos, vuale is one of types: 1=seq error(<min_heter), 2=heter snp, 3=homo snp(>min_heter). 
    Init when traverse reads in BAM file and determine when traverse SNPs in posList 
    '''

    # init a list for SNPs #
    PL = posList.PositionalList() # initialize a positional list
    p = PL.add_first(Read('Begin', 0, 0, 0, None)) # initialize the first item in PL list and store the position in varible p

    # Identify SAM/BAM file to open different pysam IO handle. #
    if os.path.splitext(input)[1] == '.sam': # SAM file
        bamfile = pysam.AlignmentFile(input, "r")
    elif os.path.splitext(input)[1] == '.bam': # BAM file
        bamfile = pysam.AlignmentFile(input, "rb") # file handle of BAM file
    else:
        raise IOError('Please choose a SAM/BAM file!')
    target = bamfile.fetch(chrom, reg_s, reg_e) # iterable method of target region read

    ##### Detect all SNP in all read #####
    for read in target:
        qname = read.query_name
        start = read.query_alignment_start 
        end = read.query_alignment_end 
        ave_sq = sum(read.query_qualities) / float(len(read.query_qualities)) #average sequencing quality
        snp = {} # dict
        for x in read.get_aligned_pairs(matches_only=False, with_seq=True): # tuple of read pos, ref pos, ref seq
            read_pos, ref_pos, ref_seq = x[0:3]
            if ref_seq is not None and ref_seq.islower(): # lower case means subsititution
                snp_pos = ref_pos
                snp_alt = read.query_sequence[read_pos]
                snp_qual = read.query_qualities[read_pos]
                snp.setdefault(snp_pos, [snp_alt, snp_qual]) # add value to key
                heter_snp.setdefault(snp_pos, 0) # add None type to heter_snp dict
        p = PL.add_after(p, Read(qname, start, end, ave_sq, snp)) # add SNPs of read to positional list
    bamfile.close() # file handle closed
    PL.delete(PL.first()) # move first "begin" item
    #print(PL.first().getNode().getElement().getQname())
    ##### Identify heterozygous SNP marker; seq error and homo SNP #####
    heter_snp = heterSnp.HeterSNP(PL, heter_snp, max_heter, min_heter)
    return 0

if __name__ == '__main__':
    import getopt, sys
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:o:c:u:w:p:d:v:s:e:m:")
    except getopt.GetoptError as err:
        usage.usage()
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    # default value
    output = None # output dirctory
    input = '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/merged5YH.best.ccs.sort.bam' # input BAM/SAM file
    temp_delete = False # if delete the temp direcory
    chrom = 'chr6' # chromesome name
    reg_s = 0 # start coordinate of the region
    reg_e = 0 # end coordinate of the region
    max_heter = 0.75 # upper heter snp cutoff, alt fre/seq depth
    min_heter = 0.24 # #lower heter snp cutoff, alt fre/seq depth
    seed_win = 300 # window size of seed region selection
    seed_cut = 0.30 # cutoff value of seed (SNP) pattern selection
    snp_coin = 0.90  # SNP coincide proportion when extending
    score_cut = 0.55 # phase 0 and 1 cutoff value of scoring
    for o, a in opts:
        if o == '-h':
            usage.usage()
            sys.exit()
        elif o == '-b':
            input = a
        elif o == '-o':
            output = a

        elif o == '-m':
            chrom = a # chromesome
        elif o == '-s':
            reg_s = int(a) # start coordinate of the region
        elif o == '-e':
            reg_e = int(a) # end coordinate of the region

        elif o == '-p':
            max_heter = a # upper heter snp cutoff, alt fre/seq depth
        elif o == '-d':
            min_heter = a # #lower heter snp cutoff, alt fre/seq depth

        elif o == '-w':
            seed_win = a # window size of seed region selection
        elif o == '-v':
            seed_cut = a # cutoff value of seed (SNP) pattern selection

        elif o == '-c':
            snp_coin = a # SNP coincide proportion when extending
        elif o == '-u':
            score_cut = a # phase 0 and 1 cutoff value of scoring
    # Run the program
    main(input, output, chrom, reg_s, reg_e, seed_win, seed_cut,  max_heter, min_heter, snp_coin, score_cut)
