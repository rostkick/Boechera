from glob import glob
import os
import pathlib
from Bio import SeqIO
g = glob('*.fasta')

pathlib.Path(os.getcwd()+'/fa7to28').mkdir(exist_ok=True)

for i in g:
    with open(i) as file:
        for record in SeqIO.parse(file, 'fasta'):
            seq_name = i.split('.')[0]

            # print(record.id)
            # print(record.seq)
            with open('fa7to28/{}.fa'.format(record.id), 'a') as outfile:
                outfile.write('>{}\n{}\n'.format(seq_name, record.seq))



















    
