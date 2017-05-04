#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' Read class'''

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

##### Data structure #####
class Read(object):
    ''' Store all usefull information of a read. (one line in BAM file.)'''
    __slots__ = '__qname', '__start', '__end', '__ave_sq', '__snp' #streamline memeory usage
    def __init__(self, qname=None, start=0, end=0, ave_sq=0, snp=None):
        self.__qname = qname # QNAME of read
        self.__start = start # start coordinate in reference genome
        self.__end = end # end coordinate in reference genome
        self.__ave_sq = ave_sq # average value of sequencing qulity
        self.__snp = snp # dict of snp. key is position, vaule is [snp_alt, snp_qual]

    def getQname(self):
        ''' Return qname of read.'''
        return self.__qname

    def getStart(self):
        ''' Return alignment start position of read.'''
        return self.__start

    def getEnd(self):
        ''' Return alignment end position of read.'''
        return self.__end

    def getAveSq(self):
        ''' Return average sequencing quality of read.'''
        return self.__ave_sq

    def getSnp(self):
        ''' Return SNPs of read.'''
        return self.__snp

    def setInfo(self, qname, start, end, ave_sq, snp):
        ''' Set information of read.'''
        self.__qname = qname
        self.__start = start
        self.__end = end
        self.__ave_sq = ave_sq
        self.__snp = snp
