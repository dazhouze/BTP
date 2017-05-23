#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''Look up sorted bam file to find --genome-size parameter by Qname (reads' ID).'''

def fetch(work_dir, bam, chrom='chr6'):
    import pysam
    import os
    if os.path.exists(work_dir) is False:
        raise ValueError('%s is not exists.' % work_dir)
    if os.path.exists(bam) is False:
        raise ValueError('%s is not exists.' % bam)

    qname_path = os.path.abspath(work_dir)
    phase_0_qname=[] # 2d list for qname
    phase_1_qname=[] # 2d list for qname
    num = 0 # phase file No.
    total_read = 0
    while os.path.exists('%s/phase_0.%d.txt' % (qname_path, num)) and os.path.exists('%s/phase_1.%d.txt' % (qname_path, num)):
        print('%s/phase_0.%d.txt' % (qname_path, num), '%s/phase_1.%d.txt' % (qname_path, num))
        phase_0_qname.append(None)
        phase_1_qname.append(None)
        qname_0 = {}
        qname_1 = {}
        with open('%s/phase_0.%d.txt' % (qname_path, num), 'r') as f0_f:
            for line in f0_f:
                 qname = line.rstrip()
                 qname_0.setdefault(qname, 1)

        with open('%s/phase_1.%d.txt' % (qname_path, num)) as f1_f:
            for line in f1_f:
                 qname = line.rstrip()
                 qname_1.setdefault(qname, 1)

        phase_0_qname[num] = qname_0
        phase_1_qname[num] = qname_1
        num += 1

    phase_0_out = [None]*num # output fastq file
    phase_1_out = [None]*num # output fastq file
    for i in range(0, num):
        phase_0_out[i] = open('%s/phase_0.%d.pos' % (qname_path, i), 'w')
        phase_1_out[i] = open('%s/phase_1.%d.pos' % (qname_path, i), 'w')
         
    bamfile = pysam.AlignmentFile(bam, "rb") # file handle of BAM file
    target = bamfile.fetch(chrom) # iterable method of target region read

    for read in target:
        qname = read.query_name
        start = read.reference_start
        end = read.reference_end
        for i in range(0, num):
            if qname in phase_0_qname[i]:
                phase_0_out[i].write('%d\t%d\n' % (start, end))
            if qname in phase_1_qname[i]:
                phase_1_out[i].write('%d\t%d\n' % (start, end))

    for i in range(0, num):
        phase_0_out[i].close()
        phase_1_out[i].close()

if __name__ == '__main__':
    '''Set default value.'''
    import sys
    work_dir = None # work_dir phased read qname folder
    bam = None # bam fastq files list
    work_dir = sys.argv[1]
    bam = sys.argv[2]
    print(' - Start look up:\n', bam)
    fetch(work_dir, bam)
    print(' - Finish fetch fastq to dirctory', work_dir)
