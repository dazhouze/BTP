#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def Evaluation(phase_0, phase_1, pos_level, read_queue, heter_snp, reg_s, reg_e, log):
    import math, os.path
    hit_p = os.path.join(log, 'hit.txt') # path of log.txt
    hit_f = open(hit_p, 'w')
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
                if snp_alt == phase_1[ind]:
                    cons_1_v += snp_qual
                    cons_1_n += 1
                if snp_alt is None: # Deletion varient in read
                    pass
                elif snp_alt != phase_0[ind] and snp_alt != phase_1[ind]:
                        conf_n += 1

        if cons_0_n == cons_1_n == conf_n == 0 and reg_s < (start+end)/2 < reg_e: # no marker
            phase_0_q.append(qname)
            phase_1_q.append(qname)
        if (cons_0_n + cons_1_n + conf_n) > 0.1*len(heter_snp) and (cons_0_n + cons_1_n + conf_n) > 2:
            if cons_0_n > cons_1_n and (cons_1_n+conf_n) < 0.1*len(heter_snp):
                phase_0_q.append(qname)
            if cons_1_n > cons_0_n and (cons_0_n+conf_n) < 0.1*len(heter_snp):
                phase_1_q.append(qname)
            hit_f.write('%d %d %d %d %d %d\n' % (cons_0_n, cons_1_n, conf_n, int(100*cons_0_v), int(100*cons_1_v), int(100*conf_v)))
    hit_f.close()
    return phase_0_q, phase_1_q
