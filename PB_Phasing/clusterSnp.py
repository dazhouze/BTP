#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def Clustering(tree, read_queue, bak_queue, heter_snp, chrom, reg_s, reg_e, tree_p, heter_p):
    '''Clustering heterozygous SNP marker by optimised binary tree algorithm.
    Use slice tree.
    Back up read queue.
    Clean homo snp tree.
    '''
    pos_level = {} # start from 1. tree level dict: k is sorted heter_snp position, v is level in tree
    level_pos = {} # start from 1. k, v reverse of pos_level
    i = 0 # index of tree level
    for k in sorted(heter_snp):
        i += 1
        pos_level.setdefault(k, i)
        level_pos.setdefault(i, k)
        level_pos[i]=k
    i = None # mem release
    del i # mem release

    #walk_len = 5 # local tree walk step length
    walk_len = 4 # local tree walk step length
    assert walk_len > 2, 'walk_len must larger than 2'
    assert walk_len < len(heter_snp), 'walk_len must small than slice tree len'
    level_s = 1
    while level_s < len(heter_snp):
        ''' Every time grow 5 level.'''
        # determine the start level and end level of each round
        level_e = level_s + walk_len - 1 # new level end
        if level_e >  len(heter_snp):
            level_e = len(heter_snp)
        # level_s go back if last round not long enough
        if level_e == len(heter_snp):
            level_s = level_e - walk_len + 1
        print(' - Level:%d-%d Read:%d Marker:%d' % (level_s, level_e, len(read_queue), len(level_pos)))
        #tree.preorder_indent(tree.root())

        tree.setdefault(level_e, 0.0) # value must be 0
        # determine heter-snp-marker pattern of each read
        cursor = read_queue.first() # read queue first Position
        while cursor != None:
            ele = cursor.getElement() # element
            start, end = ele.getStart(), ele.getEnd()
            if last_min(pos_level, start) > level_e: # already go through 5 level
                break
            # read queue back up to reduce traverse time
            if  end < level_pos[level_s]: # reads end before of the 5 level need to be delete
                next_c = read_queue.after(cursor) # cursor point to next node
                node = read_queue.delete(cursor)
                bak_queue.add_last(node)
                cursor = next_c                
                continue
            read_snp = ele.getSnp() # read heter-snp-marker dict
            ave_sq = ele.getAveSq()

            pat = pattern(start, end, read_snp, heter_snp, ave_sq, level_s, level_e, level_pos) # heter snp pattern
            #print(pat)
            for d in range(level_s, level_e+1): # exculde root level (index 0)
                x = d - level_s + 1 # index of pat array
                prev_b = pat[x-1][0] # previous snp base (0/1)
                now_b = pat[x][0] # this snp base (0/1)
                v = pat[x][1] # this snp_qual
                c = 1 # cross over
                if (prev_b==0 or prev_b==1) and (now_b==0 or now_b==1): # parent's left  node add one
                    if now_b==0:
                        tree.add_value_left(d-1, v, prev_b) # find depth-1, add_value_left means add left to depth x-1 node
                        # cross over
                        if prev_b == 0:
                            tree.add_cross_right(d-1, c, 1) # cross over
                        else:
                            tree.add_cross_right(d-1, c, 0) # cross over
                    elif now_b==1:
                        tree.add_value_right(d-1, v, prev_b)
                        # cross over
                        if prev_b == 0:
                            tree.add_cross_left(d-1, c, 1) # cross over
                        else:
                            tree.add_cross_left(d-1, c, 0) # cross over
            cursor = read_queue.after(cursor) # cursor point to next node
            # end of cursor traverse read_queue
        # clean alignment error snps and ambiguous snps
        tree.pruning(tree.root(), tree_p) # alignment error and snp error check
        #tree.homo_check(tree.root()) # homo snp check
        mis_level = tree.clean(level_s) # clean tree
        if mis_level is not None:
            heter_snp, pos_level, level_pos = level_clean(mis_level, heter_snp, pos_level, level_pos) # alignment error in heter_snp dict
            level_s = mis_level - 1
        else:
            level_s = level_e # renew level_s equall to next level_s
        # end of walk through heter_snp

    tree.preorder_indent(tree.root())
    phase_0, phase_1 = tree.linkage_result()
    #print(phase_0, phase_1)

    # move bak_queue to read_queue
    cursor = bak_queue.first()
    while cursor != None:
        next_c = bak_queue.after(cursor) # cursor point to next node
        node = bak_queue.delete(cursor)
        read_queue.add_last(node)
        cursor = next_c
    assert len(bak_queue) == 0, 'bak queue is not clean up'

    # link level pos and base
    assert len(heter_snp) == len(phase_0), 'length not equal'
    for k in range(1, len(phase_0)+1): # k is level, v is pos
        v = level_pos[k]
        phase_0[k-1] = heter_snp[v][phase_0[k-1]]
        phase_1[k-1] = heter_snp[v][phase_1[k-1]]

    with open(heter_p, 'a') as heter_f:
        heter_f.write('\n***\nHeterzygous SNP Markers Phaseing Result:\n')
        heter_f.write('Chromosome\tPosition\tPhase_0\tPhase_1\n')
        for x in range(0, len(phase_0)):
            pos = level_pos[x+1]
            heter_f.write('%s\t%d\t%s\t%s\n' % (chrom, pos, phase_0[x], phase_1[x]))
    print(phase_0, phase_1)
    return phase_0, phase_1, pos_level, read_queue, heter_snp

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

def pattern(start, end, read_snp, heter_snp, ave_sq, level_s, level_e, level_pos):
    pat = [ [3, 0] for x in range(0, (level_e-level_s+2))] # heter snp pattern
    for x in range(level_s, level_e+1):
        k = level_pos[x]
        ind = x - level_s + 1
        snp_alt, snp_qual = read_snp.get(k, ['R', ave_sq]) # set snp is ref-allele first
        if start<= k <= end:
            if snp_alt == heter_snp[k][0]: # snp allele is the maximum allele property base
                pat[ind][0] = 0
                pat[ind][1] = snp_qual
            elif snp_alt == heter_snp[k][1]: # snp allele is the sec-max allele property base
                pat[ind][0] = 1
                pat[ind][1] = snp_qual
            elif snp_alt is None: # deletion variants(None)
                pat[ind][0] = 2
                pat[ind][1] = snp_qual
            else: # other allele or 
                pat[ind][0] = 2
                pat[ind][1] = snp_qual
                assert not isinstance(snp_alt, list), 'Get Value error'
        else: # out of read region
            #pat.append(3)
            pass
    return pat

def level_clean(mis_level, heter_snp, pos_level, level_pos):
    pos = level_pos[mis_level]
    del heter_snp[pos]
    pos_level = None
    level_pos = None
    pos_level = {}
    level_pos = {}
    i = 0 # index of tree level
    for k in sorted(heter_snp):
        i += 1
        pos_level.setdefault(k, i)
        level_pos.setdefault(i, k)
    return heter_snp, pos_level, level_pos
