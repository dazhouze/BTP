#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
SNPs clustering.
Remove homo SNP and alignment error SNPS.

heter_snp list
    0-start
tree_pointer list
    0-start

pos_level dict
    1-start
level_pos dict
    1-start
pat list:
    1-start
tree:
    1-start
'''

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def clustering(tree, read_queue, bak_queue, heter_snp, chrom, reg_s, reg_e, tree_p, heter_p, rm_p):
    '''Clustering heterozygous SNP marker by optimised binary tree algorithm.
    Use slice tree.
    Back up read queue.
    Clean homo snp tree.
    '''
    print(' - Start clustering SNPs.')
    rm_f = open(rm_p, 'w')
    # start from 1. tree level dict: k is sorted heter_snp position, v is level in tree
    pos_level = {}
    # start from 1. k, v reverse of pos_level
    level_pos = {}
    i = 0 # index of tree level
    for k in sorted(heter_snp):
        i += 1
        pos_level.setdefault(k, i)
        level_pos.setdefault(i, k)
        level_pos[i] = k
    i = None # mem release
    del i # mem release

    # arrary for tree node Position to speed up pruning and clean of tree
    tree_pointer = [[None,None] for row in range(0, len(heter_snp)+1)]
    tree_pointer[0][0] = tree.root()
    tree_pointer[1][0], tree_pointer[1][1] = tree.left(tree.root()), tree.right(tree.root())
    walk_len = 3 # local tree walk step length
    assert walk_len > 2, 'walk_len must larger than 2'
    assert walk_len < len(heter_snp), 'walk_len must small than slice tree len'
    level_s = 1
    while level_s < len(heter_snp):
        ''' Every time grow walk_len level.'''
        ''' determine the start level and end level of each round.'''
        level_e = level_s + walk_len # new level end
        if level_e > len(heter_snp):
            level_e = len(heter_snp)
        # level_s go back if last round not long enough
        print('  - Level:%d-%d Read:%d Marker:%d' % \
              (level_s, level_e, len(read_queue), len(level_pos)))

        tree.setdefault(tree_pointer[level_s][0], level_e, 0)
        tree.setdefault(tree_pointer[level_s][1], level_e, 0)
        '''Determine heter-snp-marker pattern of each read.'''
        cursor = read_queue.first() # read queue first Position
        while cursor != None:
            ele = cursor.get_element() # element
            start, end = ele.get_start(), ele.get_end()
            if last_min(pos_level, start) > level_e: # already go through 5 level
                break
            '''Read queue back up to reduce traverse time.'''
            if  end < level_pos[level_s]: # reads end before of the 5 level need to be delete
                next_c = read_queue.after(cursor) # cursor point to next node
                node = read_queue.delete(cursor)
                bak_queue.add_last(node)
                cursor = next_c
                continue
            read_snp = ele.get_snp() # read heter-snp-marker dict
            # heter snp pattern
            pat = pattern(start, end, read_snp, heter_snp, level_s, level_e, level_pos)
            #print(pat)
            for d in range(level_s, level_e+1): # exculde root level (index 0)
                x = d - level_s + 1 # index of pat array
                prev_b = pat[x-1] # previous snp base (0/1)
                now_b = pat[x] # this snp base (0/1)
                v = 1 # link num
                c = 1 # cross over
                if (prev_b == 0 or prev_b == 1) and (now_b == 0 or now_b == 1):
                    '''Add value to  parent's left node.'''
                    if now_b == 0:
                        '''Find depth-1, add_value_left means add left to depth x-1 node.'''
                        tree.add_value_left(tree_pointer[level_s][0],d-1, v, prev_b)
                        tree.add_value_left(tree_pointer[level_s][1],d-1, v, prev_b)
                        '''Crossover.'''
                        if prev_b == 0:
                            tree.add_cross_right(tree_pointer[level_s][0], d-1, c, 1) # cross over
                            tree.add_cross_right(tree_pointer[level_s][1], d-1, c, 1) # cross over
                        else:
                            tree.add_cross_right(tree_pointer[level_s][0],d-1, c, 0) # cross over
                            tree.add_cross_right(tree_pointer[level_s][1],d-1, c, 0) # cross over
                    elif now_b == 1:
                        tree.add_value_right(tree_pointer[level_s][0], d-1, v, prev_b)
                        tree.add_value_right(tree_pointer[level_s][1], d-1, v, prev_b)
                        '''Crossover.'''
                        if prev_b == 0:
                            tree.add_cross_left(tree_pointer[level_s][0],d-1, c, 1) # cross over
                            tree.add_cross_left(tree_pointer[level_s][1],d-1, c, 1) # cross over
                        else:
                            tree.add_cross_left(tree_pointer[level_s][0], d-1, c, 0) # cross over
                            tree.add_cross_left(tree_pointer[level_s][1], d-1, c, 0) # cross over
            cursor = read_queue.after(cursor) # cursor point to next node
            # end of cursor traverse read_queue
        '''Alignment error and snp error check.'''
        pointers = tree.pruning(tree_pointer[level_s][0], tree_p, 0) # left subtree pruning
        for x in range(0, len(pointers)):
            tree_pointer[level_s+x+1][0] = pointers[x]
        pointers = tree.pruning(tree_pointer[level_s][1], tree_p, 1) # right subtree pruning
        for x in range(0, len(pointers)):
            tree_pointer[level_s+x+1][1] = pointers[x]
        '''Clean alignment error snps and ambiguous snps.'''
        mis_level = tree.clean(tree_pointer[level_s], chrom, level_pos, rm_f) # clean tree
        if mis_level is not None:
            '''Delete subtree at and after mis_level.'''
            tree.delete_depth(tree_pointer[level_s][0], mis_level)
            tree.delete_depth(tree_pointer[level_s][1], mis_level)
            '''Remove homo snp and alignment error snp in heter_snp dict.'''
            heter_snp, pos_level, level_pos = level_clean\
                (mis_level, heter_snp, pos_level, level_pos, heter_p)
            tree_pointer.pop(mis_level)
            level_s = mis_level - 1
        else:
            level_s = level_e # renew level_s equall to next level_s
        # end of walk through heter_snp
    #tree.preorder_indent(tree.root()) # whole MHC region: maximum recursion depth exceeded
    print(' - Finish construct phased SNPs tree. Level: %d.' % level_e)

    read_queue = move_back(bak_queue, read_queue) # move bak_queue back to read_queue

    with open(heter_p, 'w') as heter_f:
        phase_0, phase_1, phase_pos , pos_index = tree.linkage_result\
            (len(heter_snp), heter_snp, level_pos, chrom, heter_f)
    assert level_e == len(pos_index), 'SNP pos and index not equal.'
    print(' - Finish analysis phased fragment: %d.' % len(phase_pos))
    rm_f.close()

    return phase_0, phase_1, phase_pos, pos_index, read_queue, heter_snp

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

def pattern(start, end, read_snp, heter_snp, level_s, level_e, level_pos):
    ''' Return SNPs information 0, 1 euqal to index of heter_snp.
        3 is out of range.
    '''
    pat = [3]*(level_e-level_s+2) # heter snp pattern
    for x in range(level_s, level_e+1):
        k = level_pos[x]
        ind = x - level_s + 1
        snp_alt, snp_qual = read_snp.get(k, ['R', 0]) # set snp is ref-allele first
        if start <= k <= end:
            if snp_alt == heter_snp[k][0]: # snp allele is the maximum allele property base
                pat[ind] = 0
            elif snp_alt == heter_snp[k][1]: # snp allele is the sec-max allele property base
                pat[ind] = 1
            elif snp_alt is None: # deletion variants(None)
                pat[ind] = 2
            else: # other allele or
                pat[ind] = 2
                assert not isinstance(snp_alt, list), 'Get Value error'
        else: # out of read region
            pass
    return pat

def level_clean(mis_level, heter_snp, pos_level, level_pos, heter_p):
    ''' Clean the level because of homo snp and alignemnt error.
        Clean the SNP in heter_snp, pos_level and level_pos.
    '''
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

def move_back(bak_queue, read_queue):
    '''Move bak_queue to read_queue.'''
    cursor = bak_queue.first()
    while cursor != None:
        next_c = bak_queue.after(cursor) # cursor point to next node
        node = bak_queue.delete(cursor)
        read_queue.add_last(node)
        cursor = next_c
    assert len(bak_queue) == 0, 'bak queue is not clean up'
    return read_queue
