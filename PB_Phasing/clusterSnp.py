#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def Clustering(tree, read_queue, heter_snp, chrom,log):
    import time
    pos_level = {} # tree level dict: k is sorted heter_snp position, v is level in tree
    level_pos = {} # k, v reverse of pos level
    i = 0 # index of tree level
    for k in sorted(heter_snp):
        i += 1
        pos_level.setdefault(k, i)
        level_pos.setdefault(i, k)
    i = None # mem release
    del i # mem release

    # determine heter-snp-marker pattern of each read
    cursor = read_queue.first() # read queue Position
    while cursor != None:
        x = cursor.getElement()
        start, end = x.getStart(), x.getEnd()
        read_snp = x.getSnp() # read heter-snp-marker dict
        ave_sq = x.getAveSq()
        pat = [3] # heter snp pattern
        for k in sorted(heter_snp): #
            snp_alt = read_snp.get(k, 'R') # set snp is ref-allele first
            if start<= k <= end:
                if snp_alt == heter_snp[k][0]: # snp allele is the maximum allele property base
                    pat.append(0)
                elif snp_alt == heter_snp[k][1]: # snp allele is the sec-max allele property base
                    pat.append(1)
                else: # other allele or deletion variants
                    pat.append(2)
            else: # out of read region
                pat.append(3)

        prun_d = min(rightMost(pat, 3), second_large(pos_level, start))# pruning level (find right most index 3 in the left of list pat
        tree.pruning(tree.root(), prun_d) # every read try to pruing front never covered tree
        print(pat, prun_d, len(tree))
        # determine the Position of each heter-snp-marker
        # first heter-snp-marker # level 1
        for x in range(1, len(pat)): # exculde the first heter-snp level (index 1) and root level (index 0)
            ######### local tree
            tree.setdefault(x, 0) # value must be 0
            #########
            if pat[x-1]==0 or pat[x-1]==1: # parent's left  node add one
                if pat[x]==0:
                    tree.add_value_left(x-1, 1, pat[x-1]) # find depth-1, add_value_left means add left to depth x-1 node
                    # cross over
                    if pat[x-1] == 0:
                        tree.add_value_right(x-1, 1, 1) # cross over
                    else:
                        tree.add_value_right(x-1, 1, 0) # cross over
                elif pat[x]==1:
                    tree.add_value_right(x-1, 1, pat[x-1])
                    # cross over
                    if pat[x-1] == 0:
                        tree.add_value_left(x-1, 1, 1) # cross over
                    else:
                        tree.add_value_left(x-1, 1, 0) # cross over
        cursor = read_queue.after(cursor) # read queue Position

    tree.pruning(tree.root(), len(heter_snp)+1) # final pruning
    tree.preorder_indent(tree.root())
    phase_0, phase_1 = tree.linkage_result()

    for k in sorted(level_pos): # k is level, v is pos
        v = level_pos[k]
        phase_0[k-1] = heter_snp[v][phase_0[k-1]]
        phase_1[k-1] = heter_snp[v][phase_1[k-1]]
    with open(log, 'a') as log_f:
        log_f.write('\n***\nHeterzygous SNP Markers Phaseing Result:\n')
        log_f.write('Chromosome\tPosition\tPhase_0\tPhase_1\n')
        for x in range(0, len(phase_0)):
            pos = level_pos[x+1]
            log_f.write('%s\t%d\t%s\t%s\n' % (chrom, pos, phase_0[x], phase_1[x]))
    return phase_0, phase_1, pos_level

def rightMost(l, i): # pruning level
    '''Return the index of last i in left of list l'''
    prev_ind = 0
    for k in range(0, len(l)):
        v = l[k]
        if v != i:
            return prev_ind
        prev_ind = k
    return prev_ind

def second_large(pos_level, read_s): # pruning level
    '''Return the largest k's value of dict less than v:'''
    prev_l = 1
    for k in sorted(pos_level):
        v = pos_level[k]
        if k > read_s:
            return prev_l
        prev_l = v
    return prev_l
