'''
Get raw reads info through a QNAME list file and a FASTQ file
python Qname2Fastq.py qname.txt(argv1) file.fastq(argv2)
'''

__author__ = 'Zhou Ze'
__version__ = '1.0'

import re， sys

def main():
    s = set([])
    with open(sys.argv[1], 'r') as q:
        for line in q:
            line = line.rstrip()
            s.add(line)

    with open(sys.argv[2], 'r') as fq:
        while True:
            seq_id = handle.next().strip("\n")
            seq = handle.next().strip("\n")
            seq_info = handle.next().strip("\n")
            seq_qual =handle.next().strip("\n")
            if seq_id in s:
                print(seq_id)
                print(seq)
                print(seq_info)
                print(seq_qual)

if __name__ == '__main__':
    main()
