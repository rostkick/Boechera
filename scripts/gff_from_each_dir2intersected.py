from glob import glob
from os import makedirs
from shutil import copyfile
import os

g=glob('*/*sect.gff')

try:
    makedirs('globit_out')
except:
    pass

try:
    for i in g:
        copyfile(i, f'globit_out/{os.path.basename(i)}')
except:
    pass
