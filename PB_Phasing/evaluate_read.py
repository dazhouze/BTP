#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

def evaluation(phase_0, phase_1, phase_pos, pos_index, read_queue, heter_snp, reg_s, reg_e, hit_p):
    ''' Evaluat all read in the region through phased SNPs.
    phase_0, phase_1: 2D array [base, if significant]
    '''
    print(' - Start reads evaluation.')
    hit_f = open(hit_p, 'w')
    phase_0_q, phase_1_q = [None]*len(phase_0), [None]*len(phase_0) # list for qname
    phase_None_q = [] # list for None phase reads' qname
    for i in range(0, len(phase_0)): # phased region one by one
        phase_0_q_p, phase_1_q_p = [], [] # part list for qname
        for x in read_queue: # x is Read object
            qname = x.get_qname()
            start, end = x.get_start(), x.get_end()
            if start > phase_pos[i][1] or end < phase_pos[i][0]: # read out of markers region
                continue # next read
            read_snp = x.get_snp() # read heter-snp-marker dict
            cons_0_n, cons_1_n = 0, 0 # consistence number(same as heter-snp-marker)
            conf_n = 0 # conflict number(diff from heter-snp-marker)
            '''phase_0 , phase_1 heter-snp-marker hit count.'''
            for k in heter_snp: # k is pos, v is tuple of alt
                if start <= k <= end: # heter-snp-marker is in read region
                    row, col =  pos_index[k] # phase_0 phase_1 index
                    if i != row: # snp is not in this fragment
                        continue
                    v = read_snp.get(k, ['R', 0])
                    snp_alt = v[0] # if not in read_snp consider as Ref-allele (snp_alt = v[0])
                    if snp_alt is None: # Deletion varient in read
                        continue # next snp
                    
                    # only consisder the significant SNPs
                    if snp_alt == phase_0[row][col][0]: # snp alt is phase 0
                        if phase_0[row][col][1] is True: # phase 0 is significant
                            cons_0_n += 1
                    elif snp_alt == phase_1[row][col][0]: # snp alt is other
                        if phase_1[row][col][1] is True: # phase 1 is significant
                            cons_1_n += 1
                    else: # other allele
                        conf_n += 1

            hit_f.write('%d\t%d\t%d\n' % (cons_0_n, cons_1_n, conf_n))

            if cons_0_n == cons_1_n == conf_n == 0 and reg_s < (start+end)/2 < reg_e: # no marker
                phase_0_q_p.append(qname)
                phase_1_q_p.append(qname)
            # one read at lease cover 10% of or 4 heter marker
            if cons_0_n > cons_1_n: # phase 0 read
                phase_0_q_p.append(qname)
            elif cons_1_n > cons_0_n: # phase 1 read
                phase_1_q_p.append(qname)
            else:
                phase_None_q.append(qname)
            '''
            else:
                print('==', qname, cons_1_n , cons_0_n, conf_n)
                pass
            '''
        phase_0_q[i] = phase_0_q_p
        phase_1_q[i] = phase_1_q_p
    hit_f.close()

    print(' - Finish reads evaluation.')
    for x in range(0, len(phase_0)):
        print('    fragment No.%d reads number: phase_0:%d phase_1:%d' % (x, len(phase_0_q[x]), len(phase_1_q[x])))
    print('    Not phased reads:%d' % (len(phase_None_q)))

    return phase_0_q, phase_1_q, phase_None_q
