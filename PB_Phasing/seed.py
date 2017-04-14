#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'Zhou Ze'
__version__ = '0.0.1'
def Seed(PL, heter_snp, result):
    '''Set seed.
    Return a Read object.
    '''
    heter_list = sorted(heter_snp.keys()) # heter SNP position list
    reg_start = heter_list[0]  # region start
    reg_end   = heter_list[-1] # region end
    assert reg_end-reg_start > 500,'Seed select region < 500bp.'
    for x in PL: # x is Read object
        start = x.getStart()
        end = x.getEnd()
        for k, v in x.getSnp().items(): # SNP dict k:pos v:[alt, qual]
            snp_pos = k
            snp_alt = v[0]
            snp_qual = v[1]




    return result
