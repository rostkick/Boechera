import os, re
from Bio import SeqIO

os.getcwd()
os.chdir('../input_data/exonerate_out/out_ex_parse/gff_out/extracted_fasta/genome2gene/hybrids_out')

o = os.listdir()

try:
    os.mkdir('out_28')
except:
    pass

for file in o:
    for record in SeqIO.parse(file, 'fasta'):
        with open(f'out_28/{record.id}.fasta', 'a') as w:
            w.write(f'>{file}\n'
                    f'{record.seq}\n')







