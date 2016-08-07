#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Retrun a read infomation array including: qname, start, end, length, snp[]
'''
__author__ = 'Zhou Ze'
__version__ = 1.0

import pysam
import os

class Read (object):
    '''
    Retrun a read infomation array including: qname, start, end, length, snp[]
    '''
    def __init__ (self, qname, start, end, snp = None):
        self.__qname = qname
        self.__start = start
        self.__end = end
        self.__length = end - start
        if snp == None:
            self.__snp = []
        else:
            self.__snp = snp

    def getQname (self):
        return self.__qname

    def getStart (self):
        return self.__start

    def getEnd (self):
        return self.__end
    
    def getLen (self):
        return self.__length

    def getSnp (self):
        return self.__snp

def calling (ref, seq, start, cigar):
    '''Though strings compare to call SNP'''
    pos = start #initialize the pos of seq
    snp = []

    for ci in cigar:
        if ci[0] == 0: #Match & Mismatch
            bref = ref[:ci[1]] #block reference sequence
            bseq = seq[:ci[1]] #block read sequence
            ref = ref[ci[1]:] #rest of ref
            seq = seq[ci[1]:] #rest of seq
            for i in range(ci[1]):
                r = bref[i]
                s = bseq[i]
                bpos = pos + i #base coordinate on reference genome
                if r != s:
                    snp.append((bpos, r, s))
            pos = pos + ci[0]

        elif ci[0] == 1: #Insertion
            bseq = seq[:(ci[1]-1)] #block read sequence
            seq = seq[ci[1]:]

        elif ci[0] == 2: #Deletion
            bref = ref[:(ci[1]-1)] #block reference sequence
            ref = ref[ci[1]:]
            pos = pos + ci[0]

        elif ci[0] == 4: #Soft-clipping happen only in the left or right end
            bseq = seq[:(ci[1]-1)] #block read sequence
            seq = seq[ci[1]:]
            
        elif ci[0] == 5: #Hard-clipping
            pass

        else:
            pass # 3 6 7 8  and 9

    return snp

def openFile (path, start, end):
    '''open file though pysam module.'''
    readArray = []
    if os.path.splitext(path)[1] == '.bam':
        xamfile = pysam.AlignmentFile(path, 'rb')
        iterCon = xamfile.fetch("chr6", start, end)

    else:
        xamfile = pysam.AlignmentFile(path, 'r')
        iterCon = xamfile.fetch("chr6", start, end)

    for read in iterCon:
        cigar = read.cigartuples #list of cigar field
        ref = read.get_reference_sequence() #str of ref sequence
        seq = read.query_sequence #str of read sequence
        ref = ref.upper()
        seq = seq.upper()
        start = read.reference_start #left most coordinate (except left Soft-clipping)
        end = read.reference_end #end position of read
        qname = read.query_name #QNAME of read
        snp = calling(ref, seq, start, cigar)

        readArray.append(Read(qname, start, end, snp)) #put SNPs of all read in an array

    xamfile.close()
    return readArray

def main():
    openFile ('/home/zhouze/team/YH05_PacBio/YH06_YH07/6_filtered_result/output/BAM/A.sort.bam', 29910247, 29913661)

if __name__ == '__main__':
    main()


''' 
cigar:
   *M   -   BAM_CMATCH      -   0
   *I   -   BAM_CINS        -   1
   *D   -   BAM_CDEL        -   2
    N   -   BAM_CREF_SKIP   -   3
   *S   -   BAM_CSOFT_CLIP  -   4
   *H   -   BAM_CHARD_CLIP  -   5
    P   -   BAM_CPAD        -   6
    =   -   BAM_CEQUAL      -   7
    X   -   BAM_CDIFF       -   8
    NM  -   NM tag          -   9
    note: The output is a list of (operation, length) tuples, such as [(0, 30)].

block:
    read.get_block() only get the ref-coorinate of Match/Mismatch blocks

md field:
    read.get_tag('MD')
'''
