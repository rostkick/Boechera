from glob import glob
import os
import pathlib
from Bio import SeqIO
import re


pathlib.Path(os.getcwd()+'/fa1to28').mkdir(exist_ok=True)

with open(os.getcwd()+'/all_genes_from_all.fa') as file:
        for record in SeqIO.parse(file, 'fasta'):
            # file_name = (record.id).split(':')[1]
            fasta_id = re.search('[\d\w]+', record.id).group()

            file_name = re.search('(([\d\w]+ [\d\w_]+)|([\d\w]{9}:[\d]{3}))',
                                    record.description).group()

            with open('fa1to28/{}.fasta'.format(file_name), 'a') as out_file:
                out_file.write(f'>{fasta_id}\n{file_name}\n')
