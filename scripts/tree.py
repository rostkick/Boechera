'''
Tree-building script
'''

from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
from glob import glob
from os.path import basename

import os


os.chdir('../input_data/eggnog/out/mafft/trimmed/tree')

try:
    os.makedirs('plot_tree')
except:
    pass

def plot_tree(tree, output_file='out.pdf'):
    matplotlib.rc('font', size=10)
    fig = plt.figure(figsize=(10, 20), dpi=100)
    axes = fig.add_subplot(1, 1.5, 1)
    Phylo.draw(tree, axes=axes, do_show=False)

    plt.savefig('plot_tree/'+output_file, dpi=100)

g = glob('*.dnd')

for i in g:
    tree = Phylo.read(i, 'newick')
    plot_tree(tree, basename(i)+'.png')
