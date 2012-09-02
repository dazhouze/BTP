#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def evaluation(phase_0, phase_1, pos_level, read_queue, heter_snp, reg_s, reg_e, hit_p):
    ''' Evaluat all read in the region through phased SNPs.'''
    hit_f = open(hit_p, 'w')
    phase_0_q, phase_1_q = [], []
    for x in read_queue: # x is Read object
        qname = x.get_qname()
        start, end = x.get_start(), x.get_end()
        read_snp = x.get_snp() # read heter-snp-marker dict
        cons_0_n = 0 # consistence number(same as heter-snp-marker)
        cons_1_n = 0 # consistence number(same as heter-snp-marker)
        conf_n = 0 # conflict number(diff from heter-snp-marker)
        # phase_0 , phase_1 heter-snp-marker hit count
        for k in heter_snp: # k is pos, v is tuple of alt
            if start <= k <= end: # heter-snp-marker is in read region
                v = read_snp.get(k, ['R', 0])
                snp_alt = v[0] # if not in read_snp consider as Ref-allele (snp_alt = v[0])
                #/ave_sq*Q
                ind = pos_level[k] - 1 # level of tree -1, phase_0 phase_1 index
                if snp_alt == phase_0[ind]:
                    cons_0_n += 1
                if snp_alt == phase_1[ind]:
                    cons_1_n += 1
                if snp_alt is None: # Deletion varient in read
                    pass
                elif snp_alt != phase_0[ind] and snp_alt != phase_1[ind]:
                    conf_n += 1

        hit_f.write('%d %d %d\n' % (cons_0_n, cons_1_n, conf_n))
        if cons_0_n == cons_1_n == conf_n == 0 and reg_s < (start+end)/2 < reg_e: # no marker
            #phase_0_q.append(qname)
            #phase_1_q.append(qname)
            print('no marker')
        # one read at lease cover 10% of or 3 heter marker
        if (cons_0_n + cons_1_n + conf_n) > 0.1*len(heter_snp) or (cons_0_n + cons_1_n + conf_n) > 3:
            if cons_1_n > 0.9*(cons_0_n + conf_n): # phase 1 read
                phase_1_q.append(qname)
            else:
                phase_0_q.append(qname)
    hit_f.close()
    return phase_0_q, phase_1_q
