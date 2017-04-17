#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''
Python version phasing program.
More rebost.
'''

# 1.Detect all SNP in all read 
# 2.Identify heterozygous SNP marker, seq error and homo SNP 
# 3.Detect ref-allele SNP in all read 
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
    def __init__(self, qname=None, start=0, end=0, ave_sq=0, snp=None):
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

    def setInfo(self, qname, start, end, ave_sq, snp):
        self.__qname =  qname
        self.__start = start
        self.__end = end
        self.__ave_sq = ave_sq
        self.__snp = snp

def main(input, output, chrom, reg_s, reg_e, seed_win, seed_cut,  max_heter, min_heter, snp_coin, score_cut):
    '''
    input: SAM/BAM file path
    output: output directory path
    reg_s: region start coordinate
    reg_e: region end coordinate
    '''
    import os, pysam
    from PB_Phasing import posList, heterSnp, seed

    # init output directory #
    output_dir = os.path.abspath(output)
    if not os.path.exists(output):
        os.makedirs(output_dir)
    log = os.path.join(output, 'log.txt') # path of log.txt
    with open(log, 'w') as log_f:
        log_f.write('***\nTarget:\n%s:%d-%d\n' % (chrom, reg_s, reg_e))
        log_f.write('input:%s\noutput:%s\nchr:%s, start:%d, end:%d\nseed_win:%d, seed_cut:%.2f\nmax_heter:%.2f, min_heter:%.2f\nsnp_coin:%.2f, score_cut:%.2f\n' % (input, output, chrom, reg_s, reg_e, seed_win, seed_cut,  max_heter, min_heter, snp_coin, score_cut))

    # init varibls for Phasing #
    heter_snp = {} 
    ''' All possible position of SNPs.
    key is pos, vuale is one of types: 1=seq error(<min_heter), 2=heter snp, 3=homo snp(>min_heter). 
    Init when traverse reads in BAM file and determine when traverse SNPs in positional list 
    '''
    seq_depth = [0]*(reg_e-reg_s+1) # region sequencing depth accurate to base

    # init a list for SNPs #
    read_queue = posList.PositionalList() # initialize a positional list
    p = read_queue.add_first(Read('Begin', 0, 0, 0, None)) # initialize the first item in read_queue list and store the position in varible p

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
        start = read.reference_start
        end   = read.reference_end
        ave_sq = sum(read.query_qualities) / float(len(read.query_qualities)) #average sequencing quality
        snp = {} # dict
        for x in read.get_reference_positions(): # add 1 to seq_depth array
            if reg_s <= x <= reg_e: # only consider pos within target region
                seq_depth[x-reg_s] += 1
        for x in read.get_aligned_pairs(matches_only=False, with_seq=True): # tuple of read pos, ref pos, ref seq
            read_pos, ref_pos, ref_seq = x[0:3]
            if ref_seq is not None and reg_s<=ref_pos<=reg_e and ref_seq.islower(): # indels' ref_seq is None. lower case means subsititution
                snp_pos = ref_pos
                snp_alt = read.query_sequence[read_pos]
                if snp_alt is  None:
                    raise ValueError(' alt is None') 
                snp_qual = read.query_qualities[read_pos]
                snp.setdefault(snp_pos, [snp_alt, snp_qual]) # add value to key
                heter_snp.setdefault(snp_pos, 0) # add None type to heter_snp dict
        p = read_queue.add_after(p, Read(qname, start, end, ave_sq, snp)) # add SNPs of read to positional list
    bamfile.close() # file handle closed
    read_queue.delete(read_queue.first()) # move first "begin" item, by test: 
    assert read_queue.first().getNode().getElement().getQname() != 'Begin', 'Fisrt item is not removed.'

    ##### Identify heterozygous SNP marker; seq error and homo SNP (within block) #####
    heter_snp = heterSnp.HeterSNP(read_queue, heter_snp, seq_depth, reg_s, reg_e, max_heter, min_heter, log) 
    print(heter_snp)
    # seq_depth mem realse

    ##### Set seed #####
    seed_0, seed_1 = seed.Seed(read_queue, heter_snp, reg_s, reg_e, seed_win, Read(), Read(), log) # artifical seed read for 2 haplotigs
    return 0
    return 0

if __name__ == '__main__':
    import getopt, sys
    from PB_Phasing import usage
    try:
        opts, args = getopt.getopt(sys.argv[1:], "htb:o:c:u:w:p:d:v:s:e:m:")
    except getopt.GetoptError as err:
        usage.Usage()
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    # default value
    output = 'test' # output dirctory
    input = '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/merged5YH.best.ccs.sort.bam' # input BAM/SAM file
    temp_delete = False # if delete the temp direcory
    chrom = 'chr6' # chromosome name
    reg_s = 28476797 # start coordinate of the region
    reg_e = 33449354 # end coordinate of the region
    max_heter = 0.75 # upper heter snp cutoff, alt fre/seq depth
    min_heter = 0.25 # #lower heter snp cutoff, alt fre/seq depth
    seed_win = 300 # window size of seed region selection
    seed_cut = 0.30 # cutoff value of seed (SNP) pattern selection
    snp_coin = 0.90  # SNP coincide proportion when extending
    score_cut = 0.55 # phase 0 and 1 cutoff value of scoring
    for o, a in opts:
        if o == '-h':
            usage.Usage()
            sys.exit()
        elif o == '-b':
            input = a
        elif o == '-o':
            output = a

        elif o == '-m':
            chrom = a # chromosome
        elif o == '-s':
            reg_s = int(a) # start coordinate of the region
        elif o == '-e':
            reg_e = int(a) # end coordinate of the region

        elif o == '-p':
            max_heter = float(a) # upper heter snp cutoff, alt fre/seq depth
        elif o == '-d':
            min_heter = float(a) # lower heter snp cutoff, alt fre/seq depth

        elif o == '-w':
            seed_win = int(a) # window size of seed region selection
        elif o == '-v':
            seed_cut = float(a) # cutoff value of seed (SNP) pattern selection

        elif o == '-c':
            snp_coin = float(a)  # SNP coincide proportion when extending
        elif o == '-u':
            score_cut = float(a) # phase 0 and 1 cutoff value of scoring
    # Run the program
    usage.check(input, output, chrom, reg_s, reg_e, seed_win, seed_cut,  max_heter, min_heter, snp_coin, score_cut)
    main(input, output, chrom, reg_s, reg_e, seed_win, seed_cut,  max_heter, min_heter, snp_coin, score_cut)
