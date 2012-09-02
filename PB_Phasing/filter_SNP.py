#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' Filter sequencing error SNPs.'''

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

class Base(object):
    '''One position's SNP information: A, T, G, C 4 allele.'''
    __slots__ = '__a', '__c', '__g', '__t', '__d' #streamline memeory usage
    def __init__(self, a=0, c=0, g=0, t=0, d=0):
        self.__a = a
        self.__c = c
        self.__g = g
        self.__t = t
        self.__d = d # Deletion number

    def addA(self):
        ''' Add 1 to allele.'''
        self.__a += 1

    def addC(self):
        ''' Add 1 to allele.'''
        self.__c += 1

    def addG(self):
        ''' Add 1 to allele.'''
        self.__g += 1

    def addT(self):
        ''' Add 1 to allele.'''
        self.__t += 1

    def addD(self):
        ''' Add 1 to allele.'''
        self.__d += 1

    def getA(self):
        ''' Return allele count.'''
        return self.__a

    def getC(self):
        ''' Return allele count.'''
        return self.__c

    def getG(self):
        ''' Return allele count.'''
        return self.__g

    def getT(self):
        ''' Return allele count.'''
        return self.__t

    def getD(self): # Deletions number
        ''' Return allele count.'''
        return self.__d

    def count(self):
        '''Return the sum of A C G T allele number.'''
        return self.__a + self.__c + self.__g + self.__t

    def max2(self, dep):
        '''Return a tuple of the 2 maximum allele(A C G T Ref) of this position.'''
        ref = dep - self.count() # reference allele frequence
        bases = [self.__a, self.__c, self.__g, self.__t, ref] # base list
        first = bases.index(max(bases)) # largest item index
        bases[first] = -1 # "remove" the largest item
        second = bases.index(max(bases)) # get the second largest item
        base_c = ['A', 'C', 'G', 'T', 'R'] # base code
        return (base_c[first], base_c[second])

def remove(read_queue, heter_snp, seq_depth, chrom, reg_s, reg_e, max_heter, min_heter, snp_p):
    '''Retrun a dict of heter SNP marker positions.
    read_queue is the SNP positional list
    heter_snp is the candidate heter SNP position.
    '''
    snp_sum = {} # SNP info(kind and frequence) summary by position
    for x in read_queue: # x is Read object
        for k, v in x.get_snp().items(): # SNP dict k:pos v:[alt, qual]
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
            elif snp_alt is None:
                snp_sum[snp_pos].addD()
            else:
                raise ValueError('Base %s is not in ACGT.' % snp_alt)

    # determine the SNP types: 1=seq error, 2=heter SNP, 3=homo SNP, 0=unknown
    print(' - Classify deteced SNPs.')
    for x in heter_snp: # candidate heter snp positions
        ref_ap = 1 - snp_sum[x].count()/seq_depth[x-reg_s] # ref allele property
        a_ap = snp_sum[x].getA()/seq_depth[x-reg_s] # Base A allele property of alignment coverage
        c_ap = snp_sum[x].getC()/seq_depth[x-reg_s] # Base C allele property of alignment coverage
        g_ap = snp_sum[x].getG()/seq_depth[x-reg_s] # Base G allele property of alignment coverage
        t_ap = snp_sum[x].getT()/seq_depth[x-reg_s] # Base T allele property of alignment coverage
        d_ap = snp_sum[x].getD()/(seq_depth[x-reg_s]) # Deletion allele property of alignment coverage

        max2_bases = snp_sum[x].max2(seq_depth[x-reg_s])
        if seq_depth[x-reg_s] < 8:
            heter_snp[x] = 0 # discard
        elif ref_ap > max_heter: # Sequencing Error
            heter_snp[x] = 1
        elif max(a_ap, c_ap, g_ap, t_ap) > max_heter: # homo snp
            heter_snp[x] = 3
        elif min_heter < max(a_ap, c_ap, g_ap, t_ap) < max_heter: # heter snp
            if 'R' not in max2_bases: # double heter snp
                heter_snp[x] = 4 # alignment error
                print('    rm double-heter SNP:%d' % x)
            elif third([a_ap, c_ap, g_ap, t_ap, ref_ap]) > 0.05: # Alignment Error
                heter_snp[x] = 4 # alignment error
                print('    rm mis-aliged SNP:%d' % x)
            else:
                heter_snp[x] = 2

    result_heter = {} # heter snp marker result dict
    result_alig = {}  # alignment error result dict
    # homo, heter snp marker number; sequencing error number, discard snp number, alignment error number
    ho, he, se, dis, ae = 0, 0, 0, 0, 0
    with open(snp_p, 'w') as snp_f:
        snp_f.write('\n***\nHeterzygous SNP Markers and Homozygous SNPs :\n')
        snp_f.write('Base Code: 0=A 1=C 2=G 3=t 4=Ref.\n')
        snp_f.write('Tpye\tChromosome\tPosition\tSNP_1\tSNP_2\n')
        for k in sorted(heter_snp): # k is position and v is type 1:seq error 2:heter 3:homo
            if reg_s <= k <= reg_e: # only within region
                v = heter_snp[k]
                if v == 0: # discard snp
                    dis += 1
                elif v == 1: # sequecing error snp
                    se += 1
                elif v == 2: # heter snp
                    he += 1
                    max2_bases = snp_sum[k].max2(seq_depth[k-reg_s])
                    # return the tuple of 2 maximum allele 0:A 1:C 2:G 3:T 4:Ref
                    result_heter.setdefault(k, max2_bases)
                    snp_f.write('heter\t%s\t%d\t%s\t%s\n' % \
                                (chrom, k, max2_bases[0], max2_bases[1]))
                elif v == 3: # homo snp
                    ho += 1
                    ## return the tuple of 2 maximum alleles 0:A 1:C 2:G 3:T 4:Ref
                    #result_homo.setdefault(k, 1)
                    #snp_f.write('homo\t%s\t%d\n' % (chrom, k))
                elif v == 4: # alignment error
                    ae += 1
                    snp_f.write('align\t%s\t%d\n' % (chrom, k))

        snp_f.write('***\nPrimary SNP result\nHomo SNP: %d\nHeter SNP: %d\n\
                    Seq Error(not shown): %d\nDiscard SNP(<8x not shown): %d\n\
                    Alignmene Error: %d\n' % (ho, he, se, dis, ae))
    return result_heter # only return the heter snp marker and homo snp pos, seq error cost too much memory

def third(vs):
    ''' Return thrid maximum allele.'''
    values = vs
    first = values.index(max(values)) # largest item index
    values[first] = -1 # "remove" the largest item
    first = values.index(max(values)) # largest item index
    values[first] = -1 # "remove" the largest item
    first = values.index(max(values)) # largest item index
    return vs[first]
