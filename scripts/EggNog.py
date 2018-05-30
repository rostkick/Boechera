'''
Script should find matches with eggnog's ID and
make file EggNog.txt with meta-information
'''

import os, re
from Bio import SeqIO
from glob import glob
from collections import defaultdict
from pprint import pprint
import operator


os.chdir('../input_data/protein_panel')

os.listdir()

ort_files = glob('orthologs/*')

d = defaultdict(dict)
with open('oneline_protein.fasta') as file, open('prot.fa') as f, open('orthologs/arabidopsis_thaliana.orthologs') as atha_ort_file:
    all_prots = file.read()
    aof = atha_ort_file.read()
    for record in SeqIO.parse(f, 'fasta'):
        pattern1 = str(record.seq)
        pattern2 = re.search('(@)(?P<id>[\w\._]+)', re.search('(\n.*\n)({0})'.format(pattern1), all_prots).group()).group('id')
        pattern3 = re.search('(?P<EGG>[\w\.]*)\t.*{0}.*'.format(pattern2), aof).group('EGG')
        for i in ort_files:
            key = os.path.basename(i).split('.')[0]
            di = {key: ''}
            supper_key = f'{record.id}_{pattern3}_{pattern1} '
            d[supper_key].update(di)
            with open(i) as orthologs_file:
                for line in orthologs_file:
                    try:
                        pattern4 = re.search('(?:{})\t(?P<ids>.*)'.format(pattern3), line).group('ids')
                        d[supper_key][key] = pattern4
                    except:
                        pass


for i in d:
    for j in d[i]:
        if d[i][j] == '':
            d[i][j] = '.'

with open('EggNog.txt', 'w') as w:
    w.write('source_prot\teggnog\talyr\tatha\tbret\tbstr\tcrub\tchir\tesal\tquery_seq\n')
for k, v in d.items():
    splitter = k.split('_')
    source_prot, eggnog_id, seq = splitter[0], splitter[1], splitter[2]
    with open('EggNog.txt', 'a') as w:
        w.write(f'{source_prot}\t{eggnog_id}\t')
        sorted_v = sorted(v.items(), key=operator.itemgetter(0))
        for kk, vv in sorted_v:
            w.write(f'\t{vv}')
        w.write(f'\t{seq}\n')



