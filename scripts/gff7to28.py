from glob import glob
import os
import errno
from re import search
import pathlib
g = glob('*cuted.gff')
g
pathlib.Path(os.getcwd()+'/parsed7gffTO28').mkdir(exist_ok=True)

for i in g:
    with open(i) as data:
        for line in data:
            file = i.split('_')[0]
            prot = search('NP_[0-9.]+',line).group()
            line = line.split('\t')
            outline = [line[0],file]
            outline.extend(line[2:-1])
            with open('parsed7gffTO28t/{}.gff'.format(prot),'a') as out:
                out.write('\t'.join(outline)+'\n')
