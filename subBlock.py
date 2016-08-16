#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
main function of sub-block region coordinate detection
'''
import callSnp

winSize = 1000

def which_win (start, end):
    mid = (start + end)/2
    return int(mid / winSize)

def snp_density (array):
    densityDic = {}
    readDic = {}
    for n in array:
        win = which_win(n.getStart(), n.getEnd())
        snpNum = n.getSnpNum()
        if win in densityDic:
            densityDic[win] = densityDic[win] + snpNum
            readDic[win] = readDic[win] + n.getLen()
        else:
            densityDic[win] = snpNum
            readDic[win] = n.getLen()

    for k in sorted(densityDic):
        start = k * winSize + 1
        end = (k + 1) * winSize
        value = densityDic[k] / readDic[k] * winSize
        print('hs6\t%d\t%d\t%d'% (start, end, value))




def get_snp ():
    #read array more detail in callSnp.py
    readSnp = callSnp.openFile ('/home/zhouze/team/YH05_PacBio/YH06_YH07/2_bestHit/YH14567.best.sorted.bam', 28477797, 33448354)#whole HMC region
    #readSnp = callSnp.openFile ('/home/zhouze/team/YH05_PacBio/YH06_YH07/2_bestHit/YH14567.best.sorted.bam', 29910247, 29913661)#part HMC region
    #print('reads num: ', len(readSnp))
    snp_density(readSnp)

if __name__ == '__main__':
    get_snp()
