#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''
Fetch back fastq file by Qname (reads' ID
'''
def main(fqs, f1, f2, out1, out2):
    qname_1 = {}
    qname_2 = {}
    with open(f1, 'r') as f1_f:
        for line in f1_f:
             qname = line.rstrip()
             qname = '@'+qname
             qname_1.setdefault(qname, 1)

    with open(f2, 'r') as f2_f:
        for line in f2_f:
             qname = line.rstrip()
             qname = '@'+qname
             qname_2.setdefault(qname, 1)

    out1_f = open(out1, 'w')
    out2_f = open(out2, 'w')
    for x in fqs:
        l = 0 # line number (1-4)
        r1 = 0 # if fetch back for out1
        r2 = 0 # if fetch back for out2
        with open(x, 'r') as fq_f:
            print('Open file:', x)
            for line in fq_f:
                l += 1
                if l == 5:
                    l = 1
                    r1 = 0
                    r2 = 0
                if r1 == 1: # wite in out1_f
                    out1_f.write(line)
                if r2 == 1:
                    out2_f.write(line)
                if l == 1:
                    qname = line.rstrip()
                    if qname in qname_1:
                        out1_f.write('%s\n' % qname)
                        r1 = 1
                    if qname in qname_2:
                        out2_f.write('%s\n' % qname)
                        r2 = 1

if __name__ == '__main__':
    import os, sys
    f1 = sys.argv[1]
    f2 = sys.argv[2]
    print('Get Qname from:%s\t%s' % (f1, f2))
    ccs=['/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.1.ccs.fastq' , '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.3.ccs.fastq']
    for x in ccs:
        if not os.path.exists(x):
            raise ValueError('File not exists:%s' % x)
    out1 = os.path.dirname(f1)
    out1 = os.path.join(out1, 'phase_0.fq')
    out2 = os.path.dirname(f2)
    out2 = os.path.join(out2, 'phase_1.fq')
    main(ccs, f1, f2, out1, out2)
