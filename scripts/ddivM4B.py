from glob import glob
from Bio import SeqIO
from os.path import basename
from os import makedirs

g=glob('inp\\*')
fa = glob('*.fa')

try:
    makedirs('out')
except:
    pass

for record in SeqIO.parse(fa[0],'fasta'):
    all_feat= record.description.split(':')
    try:
        gene, spec = all_feat[1].split()
    except:
        gene, spec = all_feat[1].split()[0], ''

    if spec:
        for i in g:
            if spec in i and gene in i:
                records = list(SeqIO.parse(i, "fasta"))
                records.append(record)
                SeqIO.write(records,'out\\'+basename(i),"fasta")
                print(i,spec, gene)
                break
    else:
        for i in g:
            if gene in i:
                records = list(SeqIO.parse(i, "fasta"))
                records.append(record)
                SeqIO.write(records,'out\\'+basename(i),"fasta")
                print(i, gene)
                break

    
    
