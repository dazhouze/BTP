#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
main function of sub-block region coordinate detection
'''
import callSnp
def get_snp ():
    #read array more detail in callSnp.py
    readSnp = callSnp.openFile ('/home/zhouze/team/YH05_PacBio/YH06_YH07/6_filtered_result/output/BAM/A.sort.bam', 29910247, 29913661)
    print(len(readSnp))
    print(readSnp[0].getQname())

#def main():
#    openFile ('/home/zhouze/team/YH05_PacBio/YH06_YH07/6_filtered_result/output/BAM/A.sort.bam', 29910247, 29913661)

if __name__ == '__main__':
    get_snp()
