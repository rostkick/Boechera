# EVOLUTION ANALYSIS OF GENES ASSOCIATED WITH APOMIXIS IN BRASSICACEAE FAMILY
## Authors
*Rostislav Skitchenko, Mike Raiko, Vladimir Bruchin*
## Installing
You can clone whole repo by standart command:

```
git clone https://github.com/rostkick/Boechera.git
```
Also if you need to save it without commit-log you schould use:

```
git clone —depth=1 https://github.com/rostkick/Boechera.git
```

If you want to save only one file (for example without saving tree-plots) you should copy link of favorite file,
go to https://minhaskamal.github.io/DownGit/#/home, paste it, and create Download Link.

## Project description
Apomixis is an interesting object for study because of the huge potential
of its use in plants selection. This phylogenetic research of several species
of the Brassicaceae family was performed in order to study genes associated
with asexual reproduction through seeds. The main emphasis of the work was
placed on the evolution through duplications and the identifcation of paralogs
responsible for apomixis-related functions.
## Goals  
* Perform a comparative phylogenetics assay of the genomes of seven plants  
* Find the patterns between specific genes and apomixis plant-forms  
* Find orthologues genes in other representatives of the Brassicaceae family  
* Build the trees of genes of interests
## Methods
* Aligning proteins of interests to proteoms of 7 Brassicaceae species. 
* Multialigning of groups of proteins of interest.
* NJ tree building.
* MrBayes tree building.
## Directories and out-dir files description
* alignments -- contains files received from aligning orthologues groups our proteins before and after aligning
* exonerate output -- contains dirty outputs files with prot2genome-matchs-annotations
* junk trees -- contains prior simple NJ-proteins-trees
* mrbayes_stat -- contains meta-information about statistical processing
* orthologs -- containts information about orthologues groups of proteins of interests
* scripts -- contains all scripts (python, shell) used in the project
* trees -- contains final tree-plots, and tree-nexus-files received from MrBayes and FigTree
* EggNog.txt -- contains information about the relationship between orthologous groups and our proteins
* PARALOGS_OF_INTERESTS.fasta -- sequences of paralogs probably related with apomixis
## Links to open-sources genomes and proteoms
* [*Arabidopsis thaliana*](https://www.ncbi.nlm.nih.gov/genome/?term=Arabidopsis+thaliana)
* [*Arabidopsis lyrata*](https://www.ncbi.nlm.nih.gov/genome/?term=arabidopsis+lyrata)
* [*Boechera stricta*](https://www.ncbi.nlm.nih.gov/genome/?term=Boechera+stricta)
* [*Capsella rubella*](https://www.ncbi.nlm.nih.gov/genome/498)
* [*Chardamine chirsuta*](http://chi.mpipz.mpg.de/assembly.html)
* [*Eutrema salsugineum*](https://www.ncbi.nlm.nih.gov/genome/12266)
## Links to databases
* [NCBI](https://www.ncbi.nlm.nih.gov/)
* [Uniprot](http://www.uniprot.org/)
* [TAIR](https://www.arabidopsis.org/)
* [EggNog](http://eggnogdb.embl.de/#/app/home)
## Links to tools 
* [MrBayes](http://mrbayes.sourceforge.net/)
* [FigTree](http://tree.bio.ed.ac.uk/software/figtree/)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* [Mesquite](https://mesquiteproject.wikispaces.com/)
## Aсknowledgements
I\`d like to thank everybody especially Bioinformatics Institute.
