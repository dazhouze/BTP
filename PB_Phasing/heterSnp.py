#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

class Base(object):
    '''One position's SNP information: A, T, G, C 4 allele.'''
    __slots__ = '__a', '__c', '__g', '__t' #streamline memeory usage
    def __init__(self, a=0, c=0, g=0, t=0):
        self.__a = a
        self.__c = c
        self.__g = g
        self.__t = t

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
        '''Return the sum of A C G T allele number.'''
        return (self.__a + self.__c + self.__g + self.__t)

    def max2(self, dep):
        '''Return a tuple of the 2 maximum allele(A C G T Ref) of this position.'''
        ref = dep -(self.__a + self.__c + self.__g + self.__t) # reference allele frequence
        bases = [self.__a, self.__c, self.__g, self.__t, ref] # base list
        first = bases.index(max(bases)) # largest item index
        bases[first] = -1 # "remove" the largest item
        second = bases.index(max(bases)) # get the second largest item
        base_c = ['A','C','G','T','R'] # base code
        return (base_c[first], base_c[second])

def HeterSNP(read_queue, heter_snp, seq_depth, reg_s, reg_e, max_heter, min_heter, log):
    '''Retrun a dict of heter SNP marker positions.
    read_queue is the SNP positional list
    heter_snp is the candidate heter SNP position.
    '''
    snp_sum = {} # SNP info(kind and frequence) summary by position
    for x in read_queue: # x is Read object
        for k, v in x.getSnp().items(): # SNP dict k:pos v:[alt, qual]
            snp_pos = k
            snp_alt = v[0]
            snp_qual = v[1]
            snp_sum.setdefault(snp_pos, Base()) # A C G T alleles
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
        ref_ap = 1 - snp_sum[x].count()/seq_depth[x-reg_s] # ref allele property
        a_ap = snp_sum[x].getA()/seq_depth[x-reg_s] # Base A allele property
        c_ap = snp_sum[x].getC()/seq_depth[x-reg_s] # Base C allele property
        g_ap = snp_sum[x].getG()/seq_depth[x-reg_s] # Base G allele property
        t_ap = snp_sum[x].getT()/seq_depth[x-reg_s] # Base T allele property
        #print('%.2f' % max([a_ap, c_ap, g_ap, t_ap])) # for all SNPs' allele property
        if seq_depth[x-reg_s] <= 8:
            heter_snp[x] = 0 # discard
            next
        elif ref_ap > max_heter: # seq error
            heter_snp[x] = 1
            next
        elif max(a_ap, c_ap, g_ap, t_ap) > max_heter: # homo snp
            heter_snp[x] = 3
            next
        elif min_heter < max(a_ap, c_ap, g_ap, t_ap) < max_heter: # heter snp
            heter_snp[x] = 2
            next

    result_heter = {} # heter snp marker result dict
    result_homo  = {} # homo snp result dict
    ho, he, se, dis = 0, 0, 0, 0 # homo snp number, heter snp marker number, sequencing error snp number, discard snp number
    with open(log, 'a') as log_f:
        log_f.write('\n***\nHeterzygous SNP Markers and Homozygous SNPs :\n')
        log_f.write('Chromosome is same as target chromosome.\nBase Code: 0=A 1=C 2=G 3=t 4=Ref.\n')
        log_f.write('Tpye\tPosition\tSNP_1\tSNP_2\n')
        for k in sorted(heter_snp): # k is position and v is type 1:seq error 2:heter 3:homo
            if reg_s<= k <= reg_e: # only within region
                v = heter_snp[k]
                if v == 0: # discard snp
                    dis += 1
                elif v == 1: # sequecing error snp
                    se += 1
                elif v == 2: # only within region heter snp
                    he += 1
                    max2_bases = snp_sum[k].max2(seq_depth[k-reg_s])
                    result_heter.setdefault(k, max2_bases) # return the tuple of 2 maximum allele 0:A 1:C 2:G 3:T 4:Ref
                    log_f.write('heter\t%d\t%s\t%s\n' % (k, max2_bases[0], max2_bases[1]))
                elif v == 3: # homo snp
                    ho += 1
                    result_homo.setdefault(k, 1) # return the tuple of 2 maximum alleles 0:A 1:C 2:G 3:T 4:Ref
                    log_f.write('homo\t%d\n' % k)
        log_f.write('***\nHomo SNP: %d\nHeter SNP: %d\nSeq Error(not shown): %d\nDiscard SNP(<8x not shown): %d\n' % (ho,he,se,dis))
    return result_heter, result_homo # only return the heter snp marker and homo snp pos, seq error cost too much memory
