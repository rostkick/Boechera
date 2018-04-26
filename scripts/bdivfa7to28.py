
from glob import glob
import os, re
import pathlib
from Bio import SeqIO


os.getcwd()
os.chdir('../input_data/exonerate_out/out_ex_parse/gff_out/extracted_fasta/genome2gene')
os.listdir()
# os.chdir('../out_extracted_fasta')
# pathlib.Path(os.getcwd()+'/genome2gene').mkdir(exist_ok=True)

path='Bdiv.fa'

with open(path) as file:
    for record in SeqIO.parse(file, 'fasta'):
        # seq_name = i.split('.')[0]
        splitter = re.split(':', record.id)
        seq_name = re.split('\.', path)[0]
        file_name = splitter[1]
        print(file_name, seq_name)
        with open(f'bdiv/{file_name}.fa', 'a') as outfile:
            outfile.write(f'>{seq_name}\n{record.seq}\n')




















