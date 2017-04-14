#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'Zhou Ze'
__version__ = '0.0.1'


class Base(object):
    ''' A, T, G, C 4 base'''
    __slots__ = '__a', '__c', '__g', '__t' #streamline memeory usage
    def __init__(self):
        self.__a = 0
        self.__c = 0
        self.__g = 0
        self.__t = 0

    def addA(self):
        self.__a += 1
    def addC(self):
        self.__c += 1
    def addG(self):
        self.__g += 1
    def addT(self):
        self.__t += 1

    def getA(self):
        return self.__a
    def getC(self):
        return self.__c
    def getG(self):
        return self.__g
    def getT(self):
        return self.__t
    def count(self):
        return self.__a + self.__c + self.__g + self.__t

    def max2(self):
        '''Return the maxium 2 Base in A C G T and Ref-allel.'''
        first = second = 1
        bases = [self.__a, self.__c, self.__g, self.__t, self.__r]
        first = bases.index(max(bases))
        raise ValueError('must impliment by writer.')

def HeterSNP(PL, heter_snp, max_heter, min_heter):
    '''Retrun a dict of heter SNP marker positions.
    PL is the SNP positional list
    heter_snp is the candidate heter SNP position.
    '''
    seq_dep = {} # Depth of  covered region. Set when traverse reads in BAM file
    snp_sum = {} # SNP info(kind and frequence) summary by position
    for x in PL: # x is Read object
        start = x.getStart()
        end = x.getEnd()
        for i in range(start, end+1): # log the sequencing depth
            seq_dep.setdefault(i, 0)
            seq_dep[i] += 1
        for k, v in x.getSnp().items(): # SNP dict k:pos v:[alt, qual]
            snp_pos = k
            snp_alt = v[0]
            snp_qual = v[1]
            #snp_sum.setdefault(snp_pos, Base())# A C G T base
            if snp_pos not in snp_sum:
                snp_sum[snp_pos]=Base() # A C G T base
            #print(snp_alt,snp_sum[snp_pos] )
            if snp_alt == 'A':
                snp_sum[snp_pos].addA()
            elif snp_alt == 'C':
                snp_sum[snp_pos].addC()
            elif snp_alt == 'G':
                snp_sum[snp_pos].addG()
            elif snp_alt == 'T':
                snp_sum[snp_pos].addT()
            else:
                raise ValueError('Base %s is not in ACGT.' % snp_alt)

    # determine the SNP types: 1=seq error, 2=heter SNP, 3=homo SNP, 0=unknown
    for x,v in snp_sum.items(): # candidate heter snp positions
        #print(x, v.con())
        pass
    for x,v in heter_snp.items(): # candidate heter snp positions
        #print(x,v)
        #print(snp_sum[x].a,snp_sum[x].c,snp_sum[x].g,snp_sum[x].t,'=' ,seq_dep[x])

        ref_ap = 1 - snp_sum[x].count()/seq_dep[x] # ref allel proporty
        if seq_dep[x] == 1:
            heter_snp[x] = 0 # discard
        else:
            if max_heter < ref_ap: # seq error
                heter_snp[x] = 1
            else:
                # homo snp
                if snp_sum[x].getA()/seq_dep[x] > max_heter:  # A allel proporty
                    heter_snp[x] = 3
                    next
                elif snp_sum[x].getC()/seq_dep[x] > max_heter: # C allel proporty
                    heter_snp[x] = 3
                    next
                elif snp_sum[x].getG()/seq_dep[x] > max_heter: # G allel proporty
                    heter_snp[x] = 3
                    next
                elif snp_sum[x].getT()/seq_dep[x] > max_heter: # T allel proporty
                    heter_snp[x] = 3
                    next
                # heter snp
                elif min_heter <= snp_sum[x].getA()/seq_dep[x]: # A allel proporty
                    heter_snp[x] = 2
                    next
                elif min_heter <= snp_sum[x].getC()/seq_dep[x]: # C allel proporty
                    heter_snp[x] = 2
                    next
                elif min_heter <= snp_sum[x].getG()/seq_dep[x]: # G allel proporty
                    heter_snp[x] = 2
                    next
                elif min_heter <= snp_sum[x].getT()/seq_dep[x]: # T allel proporty
                    heter_snp[x] = 2
                    next
    result = {} # the result dict
    for k, v in heter_snp.items():
        if v == 2: # only return the heter snp
            result.setdefault(k, v)
    return result # only return the heter snp

