'''
filter multifasta file with hard seq_length disbalance
should been in dir with file 'spisok.txt', that contains:

file_name1--seq_length1
file_name2--seq_length2
file_name3--seq_length3
file_name4--seq_length4
        .......
in this dir locate dir "inp" with fasta files
'''




from glob import glob
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
from os import makedirs
from os.path import basename


g=glob('inp/*')
records = []

try:
    makedirs('out')
except:
    pass
with open('res.txt','w') as data:
    pass
l, h = 100, 1500
key = {}
with open('spisok.txt') as data:
    for line in data:
        if line.rstrip():
            line = line.split('--')
            key[line[0]]=int(line[1].rstrip())
for fa in g:
    rec, records = [], []
    k = 0
    minrec = []
    minval = [100000000,100000000]
    for record in SeqIO.parse(fa,'fasta'):
        f_med = key[basename(fa)]
        k += 1
        if len(str(record.seq))<=f_med+h and len(str(record.seq))>=f_med-l:
            rec.append(record)
        else:
            h_dist = abs(len(str(record.seq))-(f_med+h))
            l_dist = abs(len(str(record.seq))-(f_med-l))
            if (h_dist<minval[0] and h_dist<minval[1]) or (l_dist<minval[0] and l_dist<minval[1]):
                minrec = record
                minval = [l_dist, h_dist]


    if len(rec) == 0:
        rec.append(minrec)
        with open('res.txt','a') as data:
            data.write('\n'.join([basename(fa)+'--'+str(f_med),'interval from '+str(f_med-l)+' to '+str(f_med+h),'all '+str(k),\
                                  'fact len '+str(len(rec[0].seq)),'distanse '+str(minval[0])+'; '+str(minval[1]),'new the nearest '+str(len(rec)),'']))
    else:
        with open('res.txt','a') as data:
            data.write('\n'.join([basename(fa)+'--'+str(f_med),'interval from '+str(f_med-l)+' to '+str(f_med+h),'all '+str(k),'new '+str(len(rec)),'']))

    SeqIO.write(rec,'out/'+basename(fa), "fasta")
