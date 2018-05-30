'''
the script should using for parsing gff in bed file, because next tool will not work with gff file
:input: gff (chrN model gene x y score strand qual .... prot_name...)
:out:   bed file (chrN x y prot_name 0 strand)
'''
import re, os
from glob import glob

os.getcwd()

try:
    os.chdir('../input_data/exonerate_out/out_ex_parse')
except FileNotFoundError:
    pass

g = glob('*.gff')

try:
    os.mkdir("gff_out")
except:
    pass

l = []
for file in g:
    file_split = re.split('\.', file)[0]
    with open(file, 'r') as f, open(f'gff_out/{file_split}.bed', 'w') as w:
        for i in f:

            s = re.split('\t| ', i)
            w.write(f'{s[0]}\t{s[3]}\t{s[4]}\t{s[12]}:{s[3]}:{s[4]}\t0\t{s[6]}\n')










