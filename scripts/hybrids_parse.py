import os, re
from Bio import SeqIO


os.getcwd()
os.listdir()
os.chdir('../input_data/exonerate_out/out_ex_parse/gff_out/extracted_fasta/genome2gene/')

try:
    os.makedirs('hybrids_out')
except:
    pass

for record in SeqIO.parse('all_genes_from_all.fa', 'fasta'):
    name = re.split(':', record.id, 1)
    print(name)
    with open(f'hybrids_out/{name[0]}', 'a') as w:
        w.write(f'>{name[1]}\n'
                f'{record.seq}\n')
    # print(record.id)
    # print(record.seq)
















