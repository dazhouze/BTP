#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''Fetch back fastq file by Qname (reads' ID).'''

import sys
import os
import getopt

def fetch(work_dir, ccs):

    if os.path.exists(work_dir) is False:
        raise ValueError('%s is not exists.' % work_dir)
    for x in ccs:
        if os.path.exists(x) is False:
            raise ValueError('%s is not exists.' % x)

    qname_path = os.path.abspath(work_dir)
    phase_0_qname=[] # 2d list for qname
    phase_1_qname=[] # 2d list for qname
    num = 0 # phase file No.
    total_read = 0
    while os.path.exists('%s/phase_0.%d.txt' % (qname_path, num)) and os.path.exists('%s/phase_1.%d.txt' % (qname_path, num)):
        print('%s/phase_0.%d.txt' % (qname_path, num), '%s/phase_1.%d.txt' % (qname_path, num))
        phase_0_qname.append(None)
        phase_1_qname.append(None)
        qname_0 = {}
        qname_1 = {}
        with open('%s/phase_0.%d.txt' % (qname_path, num), 'r') as f0_f:
            for line in f0_f:
                 qname = line.rstrip()
                 qname = '@'+qname
                 qname_0.setdefault(qname, 1)

        with open('%s/phase_1.%d.txt' % (qname_path, num)) as f1_f:
            for line in f1_f:
                 qname = line.rstrip()
                 qname = '@'+qname
                 qname_1.setdefault(qname, 1)

        phase_0_qname[num] = qname_0
        phase_1_qname[num] = qname_1
        num += 1

    phase_0_out = [None]*num # output fastq file
    phase_1_out = [None]*num # output fastq file
    for i in range(0, num):
        phase_0_out[i] = open('%s/phase_0.%d.fastq' % (qname_path, i), 'w')
        phase_1_out[i] = open('%s/phase_1.%d.fastq' % (qname_path, i), 'w')
         
    for x in ccs:
        l = 0 # line number (1-4)
        f0, f1 = None, None # file handle for phase 0, phase 1
        with open(x, 'r') as fq_f:
            print('Open file:', x)
            for line in fq_f:
                l += 1
                if l == 5:
                    l, f0, f1  = 1, None, None
                if f0 is not None:
                    phase_0_out[f0].write(line)
                if f1 is not None:
                    phase_1_out[f1].write(line)
                if l == 1:
                    qname = line.rstrip()
                    for i in range(0, num):
                        if qname in phase_0_qname[i]:
                            phase_0_out[i].write('%s\n' % qname)
                            f0 = i
                        if qname in phase_1_qname[i]:
                            phase_1_out[i].write('%s\n' % qname)
                            f1 = i

    for i in range(0, num):
        phase_0_out[i].close()
        phase_1_out[i].close()

def usage():
    info = '''
Program: fetch.py 
Version: 0.2.0
Contact: Zhou Ze <zhouze@genomocis.cn>

usage: fetch.py phased/directory  ccs/fastq/files
\t-h For help information.'''
    print(info)

if __name__ == '__main__':
    '''Set default value.'''
    work_dir = None # work_dir phased read qname folder
    ccs = [] # ccs fastq files list
    work_dir = sys.argv[1]
    ccs = sys.argv[2:]
    if len(ccs) == 0: 
        ccs=['/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.1.ccs.fastq' , '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.3.ccs.fastq']
    print(' - Start fetch fastq from:\n', ccs)
    fetch(work_dir, ccs)
    print(' - Finish fetch fastq to dirctory', work_dir)
