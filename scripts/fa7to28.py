
from glob import glob
import os, re
import pathlib
from Bio import SeqIO


os.getcwd()
os.chdir('../input_data/exonerate_out/out_ex_parse/gff_out/extracted_fasta/genome2gene')
os.listdir()
# os.chdir('../out_extracted_fasta')
# pathlib.Path(os.getcwd()+'/genome2gene').mkdir(exist_ok=True)


g = glob('*.fa')
g
for i in g:
    with open(i) as file:
        for record in SeqIO.parse(file, 'fasta'):

            # seq_name = i.split('.')[0]
            splitter = re.split(':', record.id, 1)
            seq_name = re.split('\.', i)[0]
            file_name, chords = splitter[0], splitter[1]
            with open(f'genome2gene/{file_name}.fa', 'a') as outfile:
                outfile.write(f'>{seq_name}_{chords}\n{record.seq}\n')




















