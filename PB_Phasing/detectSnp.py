#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''
Python version phasing program.
More rebost.
Binary tree.
'''

def DetectSNP(chrom, reg_s, reg_e, input, read_queue):
    '''
    input: SAM/BAM file path
    output: output directory path
    reg_s: region start coordinate
    reg_e: region end coordinate
    '''
    import os, pysam
    from PB_Phasing import myRead

    # init varibls for Phasing #
    heter_snp = {} 
    ''' All possible position of SNPs.
    key is pos, vuale is one of types: 1=seq error(<min_heter), 2=heter snp, 3=homo snp(>min_heter), 4=alignment error. 
    Init when traverse reads in BAM file and determine when traverse SNPs in positional list 
    '''
    seq_depth = [0]*(reg_e-reg_s+1) # selected region sequencing depth accurate to base
    p = read_queue.add_first(myRead.Read('Begin', 0, 0, 0, None)) # initialize the first item in read_queue list and store the position in varible p

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
        p = read_queue.add_after(p, myRead.Read(qname, start, end, ave_sq, snp)) # add SNPs of read to positional list
    bamfile.close() # file handle closed
    read_queue.delete(read_queue.first()) # move first "begin" item, by test: 
    assert read_queue.first().getNode().getElement().getQname() != 'Begin', 'Fisrt item is not removed.'

    return read_queue, heter_snp, seq_depth
