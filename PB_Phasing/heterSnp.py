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
        return (self.__a + self.__c + self.__g + self.__t)

    def max2(self, dep):
        '''Return a tuple of the 2 maxium base(A C G T Ref) of position.'''
        ref = dep -(self.__a + self.__c + self.__g + self.__t) # reference allel frequence
        bases = [self.__a, self.__c, self.__g, self.__t, ref] # base list
        first = bases.index(max(bases)) # largest item index
        bases[first] = -1 # "remove" the largest item
        second = bases.index(max(bases)) # get the second largest item
        return first, second

def HeterSNP(PL, heter_snp, seq_depth, max_heter, min_heter, reg_s, reg_e):
    '''Retrun a dict of heter SNP marker positions.
    PL is the SNP positional list
    heter_snp is the candidate heter SNP position.
    '''
    snp_sum = {} # SNP info(kind and frequence) summary by position
    for x in PL: # x is Read object
        for k, v in x.getSnp().items(): # SNP dict k:pos v:[alt, qual]
            snp_pos = k
            snp_alt = v[0]
            snp_qual = v[1]
            snp_sum.setdefault(snp_pos, Base()) # A C G T base
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
    for x in heter_snp: # candidate heter snp positions
        ref_ap = 1 - snp_sum[x].count()/seq_depth[x-reg_s] # ref allel proporty
        a_ap = snp_sum[x].getA()/seq_depth[x-reg_s] # Base A allel proporty
        c_ap = snp_sum[x].getC()/seq_depth[x-reg_s] # Base C allel proporty
        g_ap = snp_sum[x].getG()/seq_depth[x-reg_s] # Base G allel proporty
        t_ap = snp_sum[x].getT()/seq_depth[x-reg_s] # Base T allel proporty
        if seq_depth[x-reg_s] <= 8:
            heter_snp[x] = 0 # discard
            next
        elif ref_ap > max_heter: # seq error
            heter_snp[x] = 1
            next
        elif max([a_ap, c_ap, g_ap, t_ap]) > max_heter: # homo snp
            heter_snp[x] = 3
            next
        elif min_heter < max([a_ap, c_ap, g_ap, t_ap]) < max_heter: # heter snp
            #print(x, snp_sum[x].getA(),snp_sum[x].getC(),snp_sum[x].getG(),snp_sum[x].getT(),'=' ,seq_depth[x-reg_s])
            #print(x, a_ap, c_ap, g_ap, t_ap, ref_ap)
            heter_snp[x] = 2
            next

    result = {} # result dict
    for k,v in heter_snp.items(): # k is position and v is type 1:seq error 2:heter 3:homo
        if v == 2 and reg_s<= k <= reg_e: # only within region heter snp
            result.setdefault(k, snp_sum[k].max2(seq_depth[k-reg_s])) # return the tuple of 2 maxium base 0:A 1:C 2:G 3:T 4:Ref
    return result # only return the heter snp

