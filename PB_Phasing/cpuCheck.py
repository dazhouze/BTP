#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Zhou Ze'
__version__ = '0.2.0'

'''
Check cpu and memory information.
'''

def info(pid):
    ''' Use module subprocess look up process infomation. '''
    import subprocess, re
    res = subprocess.getoutput('ps ux')
    lines = re.split('\n', res)
    for line in lines:
        words = re.split('\s+', line)
        if words[1]==str(pid):
            return {'user':words[0],'pid':words[1],'cpu':words[2],'mem':words[3],'vsa':words[4],'rss':words[5],'start_time':words[9]}

    return -1

def mem(pid):
    ''' Return memory in kb. '''
    m = int(info(pid)['rss'])
    return 'PID:%s\tMemory:%d(%d kb)' % (pid, m, m/1000)

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
    print('Initial:',mem(pid))
    a = [0]*100000
    print('100k item list:',mem(pid))
    a = None
    print('a=None',mem(pid))
    a = [0]*100000
    print('100k item list:',mem(pid))
    del a
    print('del a',mem(pid))

