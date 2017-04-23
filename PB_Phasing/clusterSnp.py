#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''
Clustering heterozygous SNP marker by optimised binary tree algorithm.
'''

def Clustering(tree, read_queue, bak_queue, heter_snp, chrom, reg_s, log):
    pos_level = {} # start from 1. tree level dict: k is sorted heter_snp position, v is level in tree
    level_pos = {} # start from 1. k, v reverse of pos_level
    i = 0 # index of tree level
    for k in sorted(heter_snp):
        i += 1
        pos_level.setdefault(k, i)
        level_pos.setdefault(i, k)
    i = None # mem release
    del i # mem release

    cursor = read_queue.first() # read queue first Position
    walk = 5 # local tree walk step
    for i in range(0, (len(level_pos)-1)//walk+1): # k is level, v is pos
        ''' Every time grow 5 level.'''
        level_s = i*walk + 1 # new level start
        level_e = i*walk + walk # new level end
        if level_e > len(level_pos):
            level_e = len(level_pos)
        print('To level:%d Reads num:%d' % (level_e, len(read_queue)))
        tree.setdefault(level_e, 0) # value must be 0
        # cursor go through 5 level
        # determine heter-snp-marker pattern of each read
        while cursor != None:
            ele = cursor.getElement() # element
            start, end = ele.getStart(), ele.getEnd()
            if last_min(pos_level, start) > level_e: # already go through 5 level
                break
            if end < level_pos[level_s]: # reads end before of the 5 level need to be delete
                next_c = read_queue.after(cursor) # cursor point to next node
                node = read_queue.delete(cursor)
                bak_queue.add_last(node)
                cursor = next_c                
                continue
            read_snp = ele.getSnp() # read heter-snp-marker dict
            ave_sq = ele.getAveSq()
            pat = Pattern(start, end, read_snp, heter_snp, ave_sq) # heter snp pattern
            for x in range(level_s, level_e+1): # exculde root level (index 0)
                if (pat[x-1]==0 or pat[x-1]==1) and (pat[x]==0 or pat[x]==1): # parent's left  node add one
                    '''
                    pos1 =  level_pos[x] # heter-snp position
                    pos0 =  level_pos[x-1] # heter-snp position
                    snp_alt, snp_qual = read_snp.get(k, ['R', ave_sq]) # 
                    v = snp_qual # 1/snp_qual : without/with weight
                    '''
                    pos = level_pos[x] - reg_s
                    v = 1 # 1/snp_qual : without/with weight
                    if pat[x]==0:
                        tree.add_value_left(x-1, v, pat[x-1]) # find depth-1, add_value_left means add left to depth x-1 node
                        '''
                        # cross over
                        if pat[x-1] == 0:
                            tree.add_value_right(x-1, v, 1) # cross over
                        else:
                            tree.add_value_right(x-1, v, 0) # cross over
                        '''
                    elif pat[x]==1:
                        tree.add_value_right(x-1, v, pat[x-1])
                        '''
                        # cross over
                        if pat[x-1] == 0:
                            tree.add_value_left(x-1, v, 1) # cross over
                        else:
                            tree.add_value_left(x-1, v, 0) # cross over
                        '''
            cursor = read_queue.after(cursor) # cursor point to next node

        cursor = read_queue.first()
        tree.pruning(tree.root(), level_e+1) # every walk pruing

    #tree.preorder_indent(tree.root())
    phase_0, phase_1 = tree.linkage_result()
    print(phase_0, phase_1)

    # move bak_queue to read_queue
    cursor = bak_queue.first()
    while cursor != None:
        next_c = bak_queue.after(cursor) # cursor point to next node
        node = bak_queue.delete(cursor)
        read_queue.add_last(node)
        cursor = next_c
    assert len(bak_queue) == 0, 'bak queue is not clean up'

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
    print(phase_0, phase_1)
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

def last_min(pos_level, coor): # pruning level
    '''Return the last (the largest) k's value of dict less than v:'''
    prev_l = 1
    for k in sorted(pos_level):
        v = pos_level[k]
        if k > coor:
            return prev_l
        prev_l = v
    return prev_l

def first_max(pos_level, coor): # cursor go back level
    '''Return the first (the smallest) k's value of dict larger than v:'''
    for k in sorted(pos_level):
        v = pos_level[k]
        if k > coor:
            return v

def Pattern(start, end, read_snp, heter_snp, ave_sq):
    pat = [3] # heter snp pattern
    for k in sorted(heter_snp): #
        snp_alt, snp_qual = read_snp.get(k, ['R', ave_sq]) # set snp is ref-allele first
        if start<= k <= end:
            if snp_alt == heter_snp[k][0]: # snp allele is the maximum allele property base
                pat.append(0)
            elif snp_alt == heter_snp[k][1]: # snp allele is the sec-max allele property base
                pat.append(1)
            elif snp_alt is None: # deletion variants(None)
                pat.append(2)
                pass
            else: # other allele or 
                pat.append(2)
                assert not isinstance(snp_alt, list), 'Get Value error'
        else: # out of read region
            pat.append(3)
    return pat
