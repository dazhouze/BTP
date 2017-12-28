#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Python version phasing program.
More rebost.
Binary tree.
'''

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def main(input_file, output_file, chrom, reg_s, reg_e, seq_error):
    '''
    input_file: SAM/BAM file path
    output_file: output_file directory path
    reg_s: region start coordinate
    reg_e: region end coordinate
    '''
    import os.path
    import pysam
    import time
    from PB_Phasing import positional_list, binary_tree, detect_SNP, filter_SNP, clustering_SNP, evaluate_read

    ''' init output_file directory '''
    output_file_dir = os.path.abspath(output_file)
    if not os.path.exists(output_file_dir):
        os.makedirs(output_file_dir)
    log = os.path.join(output_file, 'log') # path of log folder
    if not os.path.exists(log):
        os.makedirs(log)
    snp_p = os.path.join(log, 'prime_snp.txt') # primary result of snps
    heter_p = os.path.join(log, 'heter_snp.txt') # heter snp in binary tree
    rm_p = os.path.join(log, 'remove_snp.txt') # remove ambiguous snp in binary tree
    tree_p = os.path.join(log, 'tree_node.txt') # binary tree node linkage and crossover
    hit_p = os.path.join(log, 'read_eval.txt') # reads evaluation
    sum_p = os.path.join(log, 'summary.txt') # path of summary.txt
    sum_f = open(sum_p, 'w') # path of summary
    now_time = time.localtime()
    sum_f.write('%s-%s-%s %s:%s\n' % (now_time[0], now_time[1], now_time[2], now_time[3], now_time[4]))
    sum_f.write('***\nOptions:\ninput_file: \"%s\"\noutput_file:\"%s\"\nchromosome=%s, start=%d, end=%d\nseq_error=%.2f\n' \
                % (input_file, output_file, chrom, reg_s, reg_e, seq_error))
    print(' - Input:  \"%s\"\n - Output: \"%s\"' % (input_file, output_file))

    ''' Entry read information '''
    read_queue = positional_list.PositionalList() # initialize a positional list
    read_queue, heter_snp, seq_depth = detect_SNP.detection(chrom, reg_s, reg_e, input_file, read_queue)
    now_time = time.localtime()
    sum_f.write('%s-%s-%s %s:%s\n' % (now_time[0], now_time[1], now_time[2], now_time[3], now_time[4]))
    sum_f.write('Detected Reads number: %d\n' % len(read_queue))
    sum_f.write('Primary candidate heterozygous SNPs: %d\n' % len(heter_snp))

    ''' Identify heterozygous SNP marker; seq error and homo SNP (within block) '''
    # heter_snp dict: k is position, v is tuple for max frequency SNP and second max frequency SNP.
    heter_snp = filter_SNP.remove(read_queue, heter_snp, seq_depth, chrom, reg_s, reg_e, seq_error, snp_p) 
    now_time = time.localtime()
    sum_f.write('%s-%s-%s %s:%s\n' % (now_time[0], now_time[1], now_time[2], now_time[3], now_time[4]))
    sum_f.write('Filtered candidate heterozygous SNPs: %d\n' % len(heter_snp))

    ''' Heterozygous SNP clustering by construct binary tree. '''
    tree = binary_tree.LinkedBinaryTree() # init a heter-snp-marker tree
    tree.add_root(binary_tree.Marker(0, 'root', 0)) # root
    tree.setdefault(tree.root(), 1, 1)
    bak_queue = positional_list.PositionalList() # a back up positional list
    tree_f = open(tree_p, 'w')
    with open(tree_p, 'w') as tree_f:
        phase_0, phase_1, phase_pos, pos_index, read_queue, heter_snp = clustering_SNP.clustering\
            (tree, read_queue, bak_queue, heter_snp, chrom, reg_s, reg_e, tree_f, heter_p, rm_p)
    now_time = time.localtime()
    sum_f.write('%s-%s-%s %s:%s\n' % (now_time[0], now_time[1], now_time[2], now_time[3], now_time[4]))
    sum_f.write('Binary tree\'s heterozygous SNPs: %d\n' % len(heter_snp))

    ''' Reads evaluation. '''
    phase_0_q, phase_1_q, phase_None_q = evaluate_read.evaluation\
        (phase_0, phase_1, phase_pos, pos_index, read_queue, heter_snp, reg_s, reg_e, hit_p)
    now_time = time.localtime()
    sum_f.write('\n%s-%s-%s %s:%s' % (now_time[0], now_time[1], now_time[2], now_time[3], now_time[4]))
    sum_f.write('Result:\nPhased fragments total number: %d\n' % len(phase_0))
    n = 0 # count of phased reads
    for x in (phase_0_q, phase_1_q):
        for i in x:
            for j in i:
                n += 1
    sum_f.write('Phased reads total number: %d\n' % n)
    sum_f.write('Not phased reads total number: %d\n' % len(phase_None_q))

    ''' Reads' Qname print out. '''
    for fragment in range(0, len(phase_0_q)):
        '''fragment caused by break point(uncovered region on reference genome).'''
        out = os.path.join(output_file, 'phase_0.%d.txt' % fragment) # path of phase_0 qname
        with open(out, 'w') as out_f:
            for x in phase_0_q[fragment]:
                out_f.write('%s\n' % x)
        out = os.path.join(output_file, 'phase_1.%d.txt' % fragment) # path of phase_1 qname
        with open(out, 'w') as out_f:
            for x in phase_1_q[fragment]:
                out_f.write('%s\n' % x)
    out = os.path.join(log, 'phase_None.txt') # path of phase_None qname
    with open(out, 'w') as out_f:
        for x in phase_None_q:
            out_f.write('%s\n' % x)

    sum_f.close()
    return 0

if __name__ == '__main__': # Run the program.
    import getopt
    import sys
    from PB_Phasing import usage
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hvb:o:p:r:s:e:m:")
    except getopt.GetoptError as err:
        usage.usage()
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)

    '''Set default value.'''
    output_file = 'output' # output_file dirctory
    input_file = '' # input_file BAM/SAM file
    chrom = 'chr6' # chromosome name
    reg_s = 28476797 # start coordinate of the region
    reg_e = 33449354 # end coordinate of the region
    seq_error = 0.25 # sequencing error snp cutoff, alt_allele/seq_depth

    '''Set command line opts value, if any.'''
    for o, a in opts:
        if o == '-h':
            usage.usage()
            sys.exit()

        if o == '-v':
            usage.usage()
            usage.author()
            sys.exit()

        elif o == '-b':
            input_file = a

        elif o == '-o':
            output_file = a

        elif o == '-m':
            chrom = a # chromosome

        elif o == '-s':
            reg_s = int(a) # start coordinate of the region

        elif o == '-e':
            reg_e = int(a) # end coordinate of the region

        elif o == '-r':
            seq_error = float(a) # upper heter snp cutoff, alt fre/seq depth

    '''Run the program.'''
    usage.check(input_file, output_file, chrom, reg_s, reg_e, seq_error)
    main(input_file, output_file, chrom, reg_s, reg_e, seq_error)
