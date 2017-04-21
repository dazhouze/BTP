#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''
Check cpu and Memory information.
'''

def info(pid):
    ''' Use module subprocess look up process infomation. '''
    import subprocess, re
    res = subprocess.getoutput('ps ux')
    lines = re.split('\n', res)
    for line in lines:
        words = re.split('\s+', line)
        if words[1]==str(pid):
            return {'user':words[0],'pid':words[1],'cpu':words[2],'Mem':words[3],'vsa':words[4],'rss':words[5],'start_time':words[9]}

    return -1

def Mem(pid):
    ''' Return Memory in kb. '''
    m = float(info(pid)['rss'])
    if m < 10*1000: # 0-10k
        return 'PID:%s\tMemory:%.2f' % (pid, m)
    elif m < 10*1000000: # 11k - 10Mb (kb)
        return 'PID:%s\tMemory:%.2f kb' % (pid, m/1000)
    elif m< 1*1000000000: # 10Mb - 10Gb(Mb)
        return 'PID:%s\tMemory:%.2f Mb' % (pid, m/1000000)
    else: # >1G
        return 'PID:%s\tMemory:%.2f Gb' % (pid, m/1000000000)

'''
(0, 'USER')
(1, 'PID')
(2, '%CPU')
(3, '%MEM')
(4, 'VSZ')
(5, 'RSS')
(6, 'TTY')
(7, 'STAT')
(8, 'START')
(9, 'TIME')
(10, 'COMMAND')
'''

if __name__ == '__main__':
    import os
    pid = os.getpid()
    print('Initial:',Mem(pid))
    a = [0]*100000
    print('100k item list:',Mem(pid))
    a = None
    print('a=None',Mem(pid))
    a = [0]*100000
    print('100k item list:',Mem(pid))
    del a
    print('del a',Mem(pid))

