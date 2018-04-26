import os, re

from Bio import SeqIO
from glob import glob

os.getcwd()
os.chdir('../input_data/exonerate_out/out_ex_parse/gff_out/extracted_fasta/genome2gene')
os.listdir()
walk = os.walk('./')
l_fasta_bdiv, l_fasta_exp = [], []
for i in walk:
    print(i)
    if i[0] == './bdiv':
        for ind in i[2]:
            l_fasta_exp.append(f'{i[0]}/{ind}')
    elif i[0] == './exp_filtred':
        for ind in i[2]:
            l_fasta_bdiv.append(f'{i[0]}/{ind}')
l_fasta_bdiv
l_fasta_exp

try:
    os.makedirs('exp_bdiv')
except:
    pass

for path_to_bdiv in l_fasta_bdiv:
    patt = re.split('\.', os.path.basename(path_to_bdiv))[0]
    for path_to_exp in l_fasta_exp:
        if re.search(patt, path_to_exp) is not None:
            # print(path_to_bdiv)
            # print(path_to_exp)
            for record in SeqIO.parse(path_to_exp, 'fasta'):
                splitter = re.split('\.', os.path.basename(path_to_exp))
                with open(f'exp_bdiv/{splitter[0]}_{splitter[1]}_combine.fasta', 'a') as w:
                    w.write(f'>{record.id}\n'
                            f'{record.seq}\n')

            for record in SeqIO.parse(path_to_bdiv, 'fasta'):
                splitter = re.split('\.', os.path.basename(path_to_exp))
                with open(f'exp_bdiv/{splitter[0]}_{splitter[1]}_combine.fasta', 'a') as w:
                    w.write(f'>{record.id}\n'
                            f'{record.seq}\n')







