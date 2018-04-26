#!/usr/bin/python3

'''
Convert a FASTA alignment to Phylip format.
Dependenies: BioPython
fasta_to_phylip --input-fasta file.fasta --output-phy file.phy
'''

from Bio import AlignIO
import argparse


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input-fasta', default='/dev/stdin')
    parser.add_argument('-o', '--output-phy', default='/dev/stdout')

    return parser.parse_args()

def main():

    args = parse_args()

    with open(args.input_fasta) as handle:
        records = AlignIO.parse(handle, 'fasta')

        with open(args.output_phy, 'w') as output_handle:
            AlignIO.write(records, output_handle, 'phylip')


if __name__ == '__main__':
    main()