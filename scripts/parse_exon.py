'''
you should use this script if you want merge many exons in gene.
input: exonerate parsed output
output: fasta files of exons
'''

import os, re
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

os.getcwd()
os.chdir('../input_data/ex_out/out')

o = os.walk('./')
for i in o:
    if i[0] == './':
        files_list = i[2]

try:
    os.makedirs('parsed')
except:
    pass

for i in files_list:
    spp = i.split('.')[0]
    with open(i) as f:
        for s in f:
            if re.search('vulgar:', s) != None:
                gene_name = re.search('AT[\w\.|\-]+', s).group()
            elif re.search('\texon\t', s) != None:
                splitter = re.split('\t', s)
                splitter[1] = spp
                new_string = '\t'.join(splitter)
                with open(f'parsed/{gene_name}', 'a') as w:
                    w.write(new_string)
            elif re.search('------------', s) != None:
                try:
                    with open(f'parsed/{gene_name}', 'a') as w:
                        w.write(s)
                except:
                    pass


os.chdir('parsed')
o = os.walk('./')
for i in o:
    if i[0] == './':
        files_list = i[2]
try:
    os.makedirs('out')
except:
    pass

for i in files_list:
    new_file_name = i.split('|')[0]+'_'+i.split('|')[1]
    with open(i) as f:
        index = 0
        for s in f:
            if re.search('exon', s) != None:
                ss = s.split('\t')
                spp = s.split('\t')[1]

                with open(f'out/{spp}_{new_file_name}_{str(index)}.bed', 'a') as w:
                    w.write(f'{ss[0]}\t{ss[3]}\t{ss[4]}\t{new_file_name}_{index}_{ss[3]}:{ss[4]}\t{ss[5]}\t{ss[6]}\n')

            if re.search('------------', s) != None:
                index += 1


os.chdir('../../../')
g1 = glob('./genome_to_exonerate_and_annotaton/*.fna')
g2 = glob('./ex_out/out/parsed/out/*')

try:
    os.makedirs('extracted_fasta')
except:
    pass

super_chr_chords = defaultdict(dict)
for i in g2:
    re_pattern = os.path.basename(i)[:4]
    with open(i) as f:
        for line in f:
            spl = line.split('\t')
            chr_name, start, end, desc = spl[0], spl[1], spl[2], spl[3]
            super_chr_chords[i].update({chr_name: (start, end, desc)})

for k, v in super_chr_chords.items():

    super_exon_dict = defaultdict(list)
    for k, v in super_chr_chords.items():
        key_name = re.search('[A-Z][\w_\.\-]+', os.path.basename(k)).group()
        splitter = os.path.basename(k).split('_')
        for i in g1:
            if splitter[0][:4] in i:
                genome = i

        records = SeqIO.to_dict(SeqIO.parse(open(genome), 'fasta'))
        exon_list = []
        for locus, x_y_desc in v.items():
            extracted_seq = str(records[locus].seq)[int(x_y_desc[0]):int(x_y_desc[1])+1]

            seq_record = SeqRecord(Seq(extracted_seq), id=locus, description=f'{re_pattern}_{x_y_desc[2]}')
            exon_list.append(seq_record)

        super_exon_dict[key_name].append(exon_list)

for k, v in super_exon_dict.items():
    with open(f'extracted_fasta/{k}.fasta', 'a') as w:
        for seq in v:
            SeqIO.write(v, w, 'fasta')











