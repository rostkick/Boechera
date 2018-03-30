'''
rename NP on fitted TAIR access nubmer from access IdList
'''
from glob import glob
from os.path import basename
from Bio import Entrez
import pandas as pd
from Bio import SeqIO
from re import search
from os import makedirs


Entrez.email = "rost20151995@gmail.com"

g=glob('*.fasta')

try:
    makedirs('out')
except:
    pass


with open ('spisok.txt') as data:
    t = data.read()

NP_list=glob('NP*')
organism = 'Viridiplantae'
for i in NP_list:
    np=basename(i).split('.')[0]
    handle_search = Entrez.esearch(db='gene', term='{0}'.format(np), retmax=20)
    records_search = Entrez.read(handle_search)
    Idrec = records_search['IdList']
    handle_summary = Entrez.esummary(db='gene', id=Idrec[0])
    record_summary = Entrez.read(handle_summary)
    AT = record_summary['DocumentSummarySet']['DocumentSummary'][0]['OtherAliases'].split(',')[0]
    try:
        res = search('{0}[ :_A-z0-9]+?.fasta'.format(AT),t).group().replace(' ','_')
        with open(i) as cop,open('out/'+res,'w') as out:
            content = cop.read()
            out.write(content)
        print(np,'_',AT)
    except:
        print ('Error',np, AT)
