#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''
Fetch back fastq file by Qname (reads' ID
'''
def main(fqs, f1):
    qname_1 = {}
    with open(f1, 'r') as f1_f:
        for line in f1_f:
             qname = line.rstrip()
             qname = '@'+qname
             qname_1.setdefault(qname, 1)

    for x in fqs:
        l = 0 # line number (1-4)
        r1 = 0 # if fetch back for out1
        with open(x, 'r') as fq_f:
            print('Open file:', x)
            for line in fq_f:
                l += 1
                if l == 5:
                    l = 1
                    r1 = 0
                if r1 == 1: # wite in out1_f
                    print(line, end = '')
                if l == 1:
                    qname = line.rstrip()
                    if qname in qname_1:
                        print(line, end = '')
                        r1 = 1

if __name__ == '__main__':
    import os, sys
    f1 = sys.argv[1]
    print('Get Qname from:', f1)
    ccs=['/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.1.ccs.fastq' , '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH01/data/m150721_013320_42266_c100807592550000001823171810291583_s1_p0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH04/data/m151029_075523_42266_c100877892550000001823193803261610_s1_p0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH05/data/m151029_121436_42266_c100877892550000001823193803261611_s1_p0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH06/data/m160301_070351_42266_c100916372550000001823203004301600_s1_X0.3.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.1.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.2.ccs.fastq', '/ifs1/ST_IM/USER/zhouze/YH_MHC_PacBio/Data/CCS/YH07/data/m160301_112718_42266_c100916372550000001823203004301601_s1_X0.3.ccs.fastq']
    for x in ccs:
        if not os.path.exists(x):
            raise ValueError('File not exists:%s' % x)
    main(ccs, f1)
