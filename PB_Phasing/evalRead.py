#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def Evaluation(phase_0, phase_1, pos_level, read_queue, heter_snp, homo_snp, log):
    import math
    phase_0_q, phase_1_q = [], []
    for x in read_queue: # x is Read object
        qname = x.getQname()
        start, end = x.getStart(), x.getEnd()
        read_snp = x.getSnp() # read heter-snp-marker dict
        ave_sq = x.getAveSq()
        cons_0_n = 0 # consistence number(same as heter-snp-marker)
        cons_1_n = 0 # consistence number(same as heter-snp-marker)
        cons_0_v = 0 # 
        cons_1_v = 0 # 
        conf_n = 0 # conflict number(diff from heter-snp-marker)
        conf_v = 0 # conflict number(diff from heter-snp-marker)
        '''
        error_n = 0 # sequencing error number
        # sequencing error snp count
        for k in read_snp: # k is pos, v is alt
            if k not in heter_snp and k not in homo_snp:
                error_n += 1
        error_n /= float(end-start)
        Q = -10 * math.log10(error_n/(1-error_n)) # sequence Quality based on Solexa pipeline
        #print (Q, ave_sq)
        '''
        # phase_0 , phase_1 heter-snp-marker hit count
        for k in heter_snp: # k is pos, v is tuple of alt
            if start <= k <= end: # heter-snp-marker is in read region
                v = read_snp.get(k, ['R', ave_sq])
                snp_alt =  v[0] # if not in read_snp consider as Ref-allele (snp_alt = v[0])
                snp_qual = v[1] # snp_qual = v[1]
                #/ave_sq*Q
                ind = pos_level[k] - 1 # level of tree -1, phase_0 phase_1 index
                if snp_alt == phase_0[ind]:
                    cons_0_v += snp_qual
                    cons_0_n += 1
                elif snp_alt == phase_1[ind]:
                    cons_1_v += snp_qual
                    cons_1_n += 1
                elif snp_alt is None: # Deletion varient in read
                    pass
                else:
                    conf_n += 1 # conflict snp(maybe indels)
                    conf_v += snp_qual
                #print(qname, start, end, k, snp_alt, phase_0[ind], phase_1[ind]) # check phase hit
        '''
        v0 = 0
        if cons_0_n != 0:
            v0 = cons_0_v/cons_0_n
        v1 = 1
        if cons_1_n != 0:
            v1 = cons_1_v/cons_1_n
        vf = 0
        if conf_n != 0:
            vf = conf_v/conf_n
        '''
        if cons_0_n == 0 and cons_1_n == 0 and conf_n == 0: # no marker
            phase_0_q.append(qname)
            phase_1_q.append(qname)
        elif 0.5*cons_0_v > cons_1_v + conf_v:
            phase_0_q.append(qname)
        elif 0.5*cons_1_v > cons_0_v + conf_v:
            phase_1_q.append(qname)
        else: # 
            pass # include in neighter haplotig
        #print(int(100*cons_0_v), int(100*cons_1_v), int(100*conf_v))
        print (len(phase_0_q), len(phase_1_q), len(read_queue))
    return phase_0_q, phase_1_q
