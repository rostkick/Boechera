'''
filter multifasta file with hard seq_length disbalance
should been in dir with file 'spisok.txt', that contains:

file_name1--seq_length1
file_name2--seq_length2
file_name3--seq_length3
file_name4--seq_length4
        .......
in this dir locate dir "inp" with fasta files
UPGRADE version also add further sequence if particular species is not represent
in file. It's was made that tree always have min 7 leafes.
'''

from glob import glob
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
from os import makedirs
import os, re
from os.path import basename

os.getcwd()
os.chdir('../input_data/exonerate_out/out_ex_parse/gff_out/extracted_fasta/genome2gene/')
os.listdir()


g=glob('exp_seqs/*.fa')
g
records = []

try:
    makedirs('exp_filtred')
except:
    pass

with open('res.txt','w') as data:
    pass

l, h = 100, 1500
key = {}
with open('spisok.txt') as data:
    for line in data:
        if line.rstrip():
            k = re.split('\.', line)[0]
            v = int(re.split('--', line)[1])
            key[k] = v
            # line = line.split('--')
            # key[line[0]] = int(line[1].rstrip())

for fa in g:
    rec, records = [], []
    k = 0
    spec_all, spec_list, minrec = [], [], []
    for record in SeqIO.parse(fa, 'fasta'):
        f_med = key[re.split('\.', basename(fa))[0]]
        k += 1
        spec = record.id.split('_')[0]
        spec_all.append(spec)
        if len(str(record.seq))<=f_med+h and len(str(record.seq))>=f_med-l:
            rec.append(record)
            spec_list.append(spec)
        else:
            h_dist = abs(len(str(record.seq))-(f_med+h))
            l_dist = abs(len(str(record.seq))-(f_med-l))+200
            minrec.append([min(h_dist, l_dist), spec, record])


    spec_all = list(set(spec_all))
    spec_list = list(set(spec_list))
    spec_list = [i for i in spec_all if i not in spec_list[:]]


    minrec = sorted(minrec, key = lambda x: x[0])
    if spec_list:
        for i in minrec:
            if spec_list:
                if i[1] in spec_list:
                    rec.append(i[2])
                    spec_list.remove(i[1])
            else:
                break

        with open('res.txt', 'a') as data:
            data.write('\n'.join([basename(fa)+'--'+str(f_med),'interval from '+str(f_med-l)+' to '+str(f_med+h), 'all '+str(k),\
                                  'new the nearest '+str(len(rec)),'']))
    else:
        with open('res.txt','a') as data:
            data.write('\n'.join([basename(fa)+'--'+str(f_med),'interval from '+str(f_med-l)+' to '+str(f_med+h), 'all '+str(k), 'new '+str(len(rec)),'']))

    SeqIO.write(rec,'exp_filtred/'+basename(fa), "fasta")
##
##
