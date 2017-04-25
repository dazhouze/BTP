#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''
Python version phasing program.
More rebost.
Binary tree.
'''

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

def main(input, output, chrom, reg_s, reg_e, max_heter, min_heter):
    '''
    input: SAM/BAM file path
    output: output directory path
    reg_s: region start coordinate
    reg_e: region end coordinate
    '''
    import os, copy, pysam
    from PB_Phasing import posList, binTree, heterSnp, clusterSnp, evalRead

    ##### init output directory #####
    output_dir = os.path.abspath(output)
    if not os.path.exists(output):
        os.makedirs(output_dir)
    log = os.path.join(output, 'log.txt') # path of log.txt
    with open(log, 'w') as log_f:
        log_f.write('***\nOptions:\ninput:%s\noutput:%s\nchr:%s, start:%d, end:%d\nmax_heter:%.2f, min_heter:%.2f\n' % (input, output, chrom, reg_s, reg_e, max_heter, min_heter))

    # init varibls for Phasing #
    heter_snp = {} 
    ''' All possible position of SNPs.
    key is pos, vuale is one of types: 1=seq error(<min_heter), 2=heter snp, 3=homo snp(>min_heter), 4=alignment error. 
    Init when traverse reads in BAM file and determine when traverse SNPs in positional list 
    '''
    seq_depth = [0]*(reg_e-reg_s+1) # selected region sequencing depth accurate to base

    # init a list for SNPs #
    read_queue = posList.PositionalList() # initialize a positional list
    p = read_queue.add_first(Read('Begin', 0, 0, 0, None)) # initialize the first item in read_queue list and store the position in varible p

    # Identify SAM/BAM file to open different pysam IO handle. #
    if os.path.splitext(input)[1] == '.bam': # BAM file
        bamfile = pysam.AlignmentFile(input, "rb") # file handle of BAM file
    else:
        raise IOError('Please choose a BAM file.(sorted by position)')
    target = bamfile.fetch(chrom, reg_s, reg_e) # iterable method of target region read

    ##### Detect all SNP in all read #####
    for read in target:
        qname = read.query_name
        start = read.reference_start
        end   = read.reference_end
        ave_sq = sum(read.query_qualities) / float(len(read.query_qualities)) #average sequencing quality
        snp = {} # dict
        # count sequencing depth of each position
        for x in read.get_reference_positions(): # add 1 to seq_depth array
            if reg_s <= x <= reg_e: # only consider pos within target region
                seq_depth[x-reg_s] += 1
        ''' 
        read_pos, ref_pos, ref_seq = x[0:3]
        Soft clipping: (0, None, None)
        Deletion: (None, 29052325, 'A')
        Insertion: (88, None, None)
        Match(SNP): (116, 29052381, 'a') (117, 29052382, 'C')
        '''
        # read information entry
        for x in read.get_aligned_pairs(matches_only=False, with_seq=True): # tuple of read pos, ref pos, ref seq
            read_pos, ref_pos, ref_seq = x[0:3]
            if ref_seq is not None and reg_s<=ref_pos<=reg_e: # Match(SNP and non SNP) and Deletions
                snp_pos = ref_pos
                if read_pos is None: # Deletion variant
                    snp_alt = None
                    snp_qual = ave_sq # 
                    assert snp_qual != 0, 'qual is 0' 
                    snp.setdefault(snp_pos, [snp_alt, snp_qual]) # add deletion to read's snp dict
                if ref_seq.islower(): # lower case means subsititution(SNP)
                    snp_alt = read.query_sequence[read_pos]
                    snp_qual = read.query_qualities[read_pos]
                    assert snp_alt != None, 'alt is None' 
                    snp.setdefault(snp_pos, [snp_alt, snp_qual]) # add value to key
                    heter_snp.setdefault(snp_pos, 0) # add None type to heter_snp dict
        p = read_queue.add_after(p, Read(qname, start, end, ave_sq, snp)) # add SNPs of read to positional list
    bamfile.close() # file handle closed
    read_queue.delete(read_queue.first()) # move first "begin" item, by test: 
    assert read_queue.first().getNode().getElement().getQname() != 'Begin', 'Fisrt item is not removed.'
    del target, read ,qname, start, end, snp, ave_sq, snp_pos, snp_alt, snp_qual

    ##### Identify heterozygous SNP marker; seq error and homo SNP (within block) #####
    heter_snp, homo_snp = heterSnp.HeterSNP(read_queue, heter_snp, seq_depth, chrom, reg_s, reg_e, max_heter, min_heter, log) # heter_snp dict: k is position, v is tuple for max frequency SNP and second max frequency SNP. homo_snp dict: k is position, v is 1
    #seq_depth = None # mem release
    #del seq_depth

    ##### Heterozygous SNP clustering by construct binary tree. #####
    tree = binTree.LinkedBinaryTree() # init a heter-snp-marker tree
    tree.add_root(binTree.Marker(0, 'root'))
    tree.setdefault(1,1)
    bak_queue = posList.PositionalList() # a back up positional list
    phase_0, phase_1, pos_level, read_queue, heter_snp = clusterSnp.Clustering(tree, read_queue, bak_queue, heter_snp, chrom, reg_s, log)
    bak_queue = None
    del bak_queue

    ##### Reads phasing. #####
    phase_0_q, phase_1_q = evalRead.Evaluation(phase_0, phase_1, pos_level, read_queue, heter_snp, output_dir)

    ##### Reads' Qname print out. #####
    out = os.path.join(output, 'phase_0.txt') # path of log.txt
    with open(out, 'w') as out_f:
        for x in phase_0_q:
            out_f.write('%s\n' % x)
    out = os.path.join(output, 'phase_1.txt') # path of log.txt
    with open(out, 'w') as out_f:
        for x in phase_1_q:
            out_f.write('%s\n' % x)

    return 0

if __name__ == '__main__':
    ''' Run the program. '''
    import getopt, sys
    from PB_Phasing import usage
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:o:p:d:s:e:m:")
    except getopt.GetoptError as err:
        usage.Usage()
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)

    # set default value
    output = 'test' # output dirctory
    input = '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/merged5YH.best.ccs.sort.bam' # input BAM/SAM file
    chrom = 'chr6' # chromosome name
    reg_s = 28476797 # start coordinate of the region
    reg_e = 33449354 # end coordinate of the region
    max_heter = 0.75 # upper heter snp cutoff, alt fre/seq depth
    min_heter = 0.25 # #lower heter snp cutoff, alt fre/seq depth

    # set command line opts value, if any
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

    # Run the program
    usage.Check(input, output, chrom, reg_s, reg_e, max_heter, min_heter)
    main(input, output, chrom, reg_s, reg_e, max_heter, min_heter)
