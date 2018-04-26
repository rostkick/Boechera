from Bio.Seq import Seq
import os
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_rna, generic_protein

os.getcwd()
os.chdir('../input_data/exonerate_out/out_ex_parse/gff_out/extracted_fasta/genome2gene/exp_bdiv/test')
os.listdir()


with open('AT1G73590_fa_combine.fasta') as f, open('0IJ59_extended_members.txt', 'a') as w:
    for record in SeqIO.parse(f, "fasta"):
        name_orf1, name_orf2, name_orf3 = f'{record.id}_orf1', f'{record.id}_orf2', f'{record.id}_orf3'
        orf1 = Seq(str(record.seq), generic_dna).transcribe().translate()
        orf2 = Seq(str(record.seq)[1:], generic_dna).transcribe().translate()
        orf3 = Seq(str(record.seq)[2:], generic_dna).transcribe().translate()
        print(name_orf1)
        print(orf1)
        print(name_orf2)
        print(orf2)
        print(name_orf3)
        print(orf3)
        w.write(f'\n>{name_orf1}\n{orf1}\n>{name_orf2}\n{orf2}\n>{name_orf3}\n{orf3}\n')






