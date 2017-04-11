#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
    return 'PID:%s\tMemory:%d kb' %(pid, int(info(pid)['rss'])/1000)

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
    print(mem(pid))
