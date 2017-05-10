#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Python version phasing program.
More rebost.
Binary tree.
'''

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def main(input, output, chrom, reg_s, reg_e, max_heter, min_heter):
    '''
    input: SAM/BAM file path
    output: output directory path
    reg_s: region start coordinate
    reg_e: region end coordinate
    '''
    import os.path
    import pysam
    from PB_Phasing import positional_list, binary_tree, detect_SNP, rm_error_SNP, clustering_SNP, evaluate_read

    ##### init output directory #####
    output_dir = os.path.abspath(output)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    log = os.path.join(output, 'log') # path of log folder
    if not os.path.exists(log):
        os.makedirs(log)
    snp_p = os.path.join(log, 'snp.txt') # path of snp.txt
    tree_p = os.path.join(log, 'tree.txt') # path of snp phasing result tree.txt
    heter_p = os.path.join(log, 'heter_snp.txt') # path of phased heter snp.txt
    hit_p = os.path.join(log, 'hit.txt') # path of log.txt
    sum_p = os.path.join(log, 'summary.txt') # path of summary.txt
    with open(sum_p, 'w') as sum_f:
        sum_f.write('***\nOptions:\ninput:%s\noutput:%s\nchr:%s, start:%d, end:%d\n\
                    max_heter:%.2f, min_heter:%.2f\n' % (input, output, chrom, reg_s,\
                                                         reg_e, max_heter, min_heter))

    ##### Entry read information #####
    read_queue = positional_list.PositionalList() # initialize a positional list
    read_queue, heter_snp, seq_depth = detect_SNP.detection(chrom, reg_s, reg_e, input, read_queue)

    ##### Identify heterozygous SNP marker; seq error and homo SNP (within block) #####
    heter_snp = rm_error_SNP.remove(read_queue, heter_snp, seq_depth, \
                                  chrom, reg_s, reg_e, max_heter, min_heter, snp_p) # heter_snp dict: k is position, v is tuple for max frequency SNP and second max frequency SNP.
    print(heter_snp)

    ##### Heterozygous SNP clustering by construct binary tree. #####
    tree = binary_tree.LinkedBinaryTree() # init a heter-snp-marker tree
    tree.add_root(binary_tree.Marker(0, 'root', 0)) # root
    tree.setdefault(1, 1)
    bak_queue = positional_list.PositionalList() # a back up positional list
    phase_0, phase_1, pos_level, read_queue, heter_snp = clustering_SNP.clustering(tree, read_queue, bak_queue, heter_snp, chrom, reg_s, reg_e, tree_p, heter_p)
    bak_queue = None
    del bak_queue

    ##### Reads phasing. #####
    phase_0_q, phase_1_q = evaluate_read.evaluation(phase_0, phase_1, pos_level, read_queue, heter_snp, reg_s, reg_e, hit_p)

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

if __name__ == '__main__': # Run the program.
    import getopt
    import sys
    from PB_Phasing import usage
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hvb:o:p:d:s:e:m:")
    except getopt.GetoptError as err:
        usage.usage()
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)

    # set default value
    output = 'HLA_test' # output dirctory
    input = '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/merged5YH.best.ccs.sort.bam' # input BAM/SAM file
    chrom = 'chr6' # chromosome name
    reg_s = 28476797 # start coordinate of the region
    reg_e = 33449354 # end coordinate of the region
    max_heter = 0.75 # upper heter snp cutoff, alt fre/seq depth
    min_heter = 0.25 # #lower heter snp cutoff, alt fre/seq depth

    # set command line opts value, if any
    for o, a in opts:
        if o == '-h':
            usage.usage()
            sys.exit()

        if o == '-v':
            usage.usage()
            usage.author()
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
    usage.check(input, output, chrom, reg_s, reg_e, max_heter, min_heter)
    main(input, output, chrom, reg_s, reg_e, max_heter, min_heter)
