#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

''' Consider that the heter snp may can not be link together.
I have to use seed and extension strategy to illustrate the heter-information.
'''

def Seed(read_queue, heter_snp, reg_s, reg_e, seed_win, result_0, result_1, log):
    '''Set seed and return 2 aritifical Read object with heter-heter-snp-marker.
    The heter-snp-marker distance: 200bp possibly cover 2 heter-heter-snp-marker
    Min.    1st Qu.    Median    Mean    3rd Qu.    Max.   
    1.0     10.0       39.0      261.6   172.0      27460.0 
    '''
    assert reg_e - reg_s > 1500,'Seed select region < 1500bp.'
    
    # set seed region with good sequencing depth and most heter-snp-marker
    # Firstly sequencing covreage:
    snp_cover = {} # heter-snp-marker coverage dict (heter-snp-marker covered by how many reads)
    for x in read_queue: # x is Read object
        start = x.getStart()
        end = x.getEnd()
        for x in range(start, end+1):
            if x in heter_snp: # only count position with heter-snp-marker
                snp_cover.setdefault(x, 0)
                snp_cover[x] += 1
    max_k, max_v = MaxDictValue(snp_cover) # max_k is the pos, max_v is the max heter-snp-marker coverage
    #print(max_k, max_v,snp_cover) # max-coverage heter-snp-marker (pos and dep)

    # find the template read: longest read which covered the max-coverage heter-snp-marker
    p = read_queue.first()
    long_l = 0 # longest read length
    long_p = None # longest read Position
    while p:
        read = p.getNode().getElement()
        r_s = read.getStart() # read start
        r_e = read.getEnd()   # read end
        r_l = r_e - r_s       # read length
        if r_s < max_k <= r_e and r_l > long_l:
            long_l = r_l
            long_p = p
        p = read_queue.after(p)
    template_read = long_p.getNode().getElement() # the template read: longest read which covered the max-coverage heter-snp-marker
    
    # construct longest artifical seed from heter-snp-marker
    # determine the start and end of seed
    tr_s = template_read.getStart() # template_read start
    tr_e = template_read.getEnd() # template_read end
    tr_snp = template_read.getSnp() # template_read heter-snp-marker dict
    seed_s = tr_s # seed start position
    seed_e = tr_e # seed end position
    for x in range(tr_s, tr_e+1):
        if x in heter_snp and snp_cover[x] < 0.55*max_v: # snp_cover: heter-snp-marker coverage max_v: max-coverage heter-snp-marker covreage(dep)
            '''heter-snp-marker coverage less than 60% max-coverage heter-snp-marker.'''
            if x < max_k: # max_k: position of max-coverage heter-snp-marker
                '''shrink left side of max-coverage heter-snp-marker.'''
                seed_s = x+1
            else:
                '''shrink right side of max-coverage heter-snp-marker.'''
                seed_e = x-1
    # set seed region heter-snp-marker
    seed_snp = {} # seed region heter_snp
    for x in range(seed_s, seed_e+1):
        if x in heter_snp:
            seed_snp.setdefault(x, heter_snp[x])
    print(tr_s, seed_s,seed_e, tr_e)

    # determine the heter-snp-marker of seed
    seed_pat = [] # seed pattern
    for x in read_queue: # x is Read object
        start = x.getStart()
        end = x.getEnd()
        read_snp = x.getSnp() # read heter-snp-marker dict
        if start<=(seed_s+seed_e)/2<=end or seed_s <=(start+end)/2<=seed_e : # read overlap with seed
            test = read_snp
            pat = []
            for k in sorted(seed_snp): #
                if start <= k <= end: # in read region
                    alt = read_snp.get(k, 'R') # set snp is ref-allele first
                    if alt == seed_snp[k][0]: # snp allele is the maximum allele property base
                        pat.append(0)
                    elif alt == seed_snp[k][1]: # snp allele is the sec-max allele property base
                        pat.append(1)
                    else: # other allele or deletion variants
                        pat.append(2)
                else: # out of read region
                    pat.append(3)
            print(pat, test==read_snp)
    return result_0, result_1
        #result_0.setInfo("seed_0", seed_s, seed_e, , seed_snp_0)
        #result_1.setInfo("seed_1", seed_s, seed_e, , seed_snp_1)

    with open(log, 'a') as log_f:
        log_f.write('\n***\nSeed i\Information:\n')
    return result_0, result_1

def MaxDictValue(dt):
    max_v = 0 # max value
    max_k = 0 # max value's key
    for k,v in dt.items():
        if v > max_v:
            max_v = v
            max_k = k
    return max_k, max_v
