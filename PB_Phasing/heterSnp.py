#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'Zhou Ze'
__version__ = '0.0.1'


class Base(object):
    ''' A, T, G, C 4 base'''
    __slots__ = 'a', 'c', 'g', 't', 'r' #streamline memeory usage
    def __init__(self):
        self.a = 0
        self.c = 0
        self.g = 0
        self.t = 0
        self.r = 0 # reference allel

    def addA(self):
        self.a += 1
    def addC(self):
        self.c += 1
    def addG(self):
        self.g += 1
    def addT(self):
        self.t += 1
    def addR(self):
        self.r += 1

def HeterSNP(PL, heter_snp, max_heter, min_heter):
    '''Retrun a dict of heter SNP marker positions.
    PL is the SNP positional list
    heter_snp is the candidate heter SNP position.
    '''
    seq_dep = {} # Depth of  covered region. Set when traverse reads in BAM file
    snp_sum = {} # SNP summary by position
    for x in PL:
        start = x.getStart()
        end = x.getEnd()
        for i in range(start, end+1): # log the sequencing depth
            seq_dep.setdefault(i, 0)
            seq_dep[i] += 1
        for k, v in x.getSnp().items(): # SNP dict k:pos v:[alt, qual]
            snp_pos = k
            snp_alt = v[0]
            snp_qual = v[1]
            snp_sum.setdefault(snp_pos, Base()) # A C G T base
            if snp_alt == 'A':
                snp_sum[snp_pos].addA()
            elif snp_alt = 'C':
                snp_sum[snp_pos].addC()
            elif snp_alt = 'G':
                snp_sum[snp_pos].addG()
            elif snp_alt = 'T':
                snp_sum[snp_pos].addT()
            else:
                raise ValueError('Base %s is not in ACGT.' % snp_alt)
    # determine the SNP types: 1=seq error, 2=heter SNP, 3=homo SNP, 0=unknown
    for x,v in heter_snp: # candidate heter snp positions
        rf = 1 - (snp_sum[x].a+snp_sum[x].c+snp_sum[x].g+snp_sum[x].t)/seq_dep[x] # ref allel proporty
        if seq_dep[x] == 1:
            v = 0 # discard
        else:
            if max_heter < rf: # seq error
                v = 1
            else:
                if snp_sum[x].a/seq_dep[x] > max_heter:
                    v = 3
                    next
                elif snp_sum[x].c/seq_dep[x] > max_heter:
                    v = 3
                    next
                elif snp_sum[x].g/seq_dep[x] > max_heter:
                    v = 3
                    next
                elif snp_sum[x].t/seq_dep[x] > max_heter:
                    v = 3
                    next
                elif min_heter <= snp_sum[x].a/seq_dep[x]: # A allel proporty
                    v = 2
                    next
                elif min_heter <= snp_sum[x].c/seq_dep[x]: # C allel proporty
                    v = 2
                    next
                elif min_heter <= snp_sum[x].g/seq_dep[x]: # G allel proporty
                    v = 2
                    next
                elif min_heter <= snp_sum[x].t/seq_dep[x]: # T allel proporty
                    v = 2
                    next
    return heter_snp

