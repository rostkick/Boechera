#!/usr/bin/R

library(ape)
library(seqinr)
library(argparse)


parser <- ArgumentParser(description='Process some integers')
parser$add_argument("-i", "--input", type="character", help='paste here your multifasta file')
parser$add_argument("-o", "--outgroup", type="character", help='outgroup name', default='Esalsugineum_173_v1_31')
parser$add_argument("-t", "--title", type="character", help='title for plot', default='Title')
parser$add_argument("-ty", "--type", type="character", help='type of tree', default='rooted')
parser$add_argument("-m", "--molecula", type="character", help='type of biopolymer DNA/protein', default='DNA')
parser$add_argument("-n", "--minpcnongap", help='minpcnogap penalty', default=30)
parser$add_argument("-g", "--minpcid", help='minpcid penalty', default=30)
args <- parser$parse_args()

# setwd('/home/rskick/IB/Projects/Boechera/input_data/prep_to_align/out/out/align/oneline/')
dir.create('out')

# args=list()
# args$input = 'AT1G02580.1_MEA.fa_mafft.fasta'
# args$title = 'Title'
# args$outgroup = 'Esalsugineum'
# args$type = 'rooted'
# args$molecula = 'DNA'
# args$minpcnongap = 30
# args$minpcid = 30



dirty_alignment <- read.alignment(file = args$input, format = 'fasta')

printMultipleAlignment <- function(alignment, chunksize=60){
    # this function requires the Biostrings package
    require("Biostrings")
    # find the number of sequences in the alignment
    numseqs <- alignment$nb
    # find the length of the alignment
    alignmentlen <- nchar(alignment$seq[[1]])
    starts <- seq(1, alignmentlen, by=chunksize)
    n <- length(starts)
    # get the alignment for each of the sequences:
    aln <- vector()
    lettersprinted <- vector()
    for (j in 1:numseqs)
    {
      alignmentj <- alignment$seq[[j]]
      aln[j] <- alignmentj
      lettersprinted[j] <- 0
    }
    # print out the alignment in blocks of 'chunksize' columns:
    for (i in 1:n) { # for each of n chunks
      for (j in 1:numseqs)
      {
        alnj <- aln[j]
        chunkseqjaln <- substring(alnj, starts[i], starts[i]+chunksize-1)
        chunkseqjaln <- toupper(chunkseqjaln)
        # Find out how many gaps there are in chunkseqjaln:
        gapsj <- countPattern("-",chunkseqjaln) # countPattern() is from Biostrings package
        # Calculate how many residues of the first sequence we have printed so far in the alignment:
        lettersprinted[j] <- lettersprinted[j] + chunksize - gapsj
        print(paste(chunkseqjaln,lettersprinted[j]))
      }
      print(paste(' '))
    }
  }
cleanAlignment <- function(alignment, minpcnongap, minpcid){
    # make a copy of the alignment to store the new alignment in:
    newalignment <- alignment
    # find the number of sequences in the alignment
    numseqs <- alignment$nb
    # empty the alignment in "newalignment")
    for (j in 1:numseqs) { newalignment$seq[[j]] <- "" }
    # find the length of the alignment
    alignmentlen <- nchar(alignment$seq[[1]])
    # look at each column of the alignment in turn:
    for (i in 1:alignmentlen)
    {
      # see what percent of the letters in this column are non-gaps:
      nongap <- 0
      for (j in 1:numseqs)
      {
        seqj <- alignment$seq[[j]]
        letterij <- substr(seqj,i,i)
        if (letterij != "-") { nongap <- nongap + 1}
      }
      pcnongap <- (nongap*100)/numseqs
      # Only consider this column if at least minpcnongap % of the letters are not gaps:
      if (pcnongap >= minpcnongap)
      {
        # see what percent of the pairs of letters in this column are identical:
        numpairs <- 0; numid <- 0
        # find the letters in all of the sequences in this column:
        for (j in 1:(numseqs-1))
        {
          seqj <- alignment$seq[[j]]
          letterij <- substr(seqj,i,i)
          for (k in (j+1):numseqs)
          {
            seqk <- alignment$seq[[k]]
            letterkj <- substr(seqk,i,i)
            if (letterij != "-" && letterkj != "-")
            {
              numpairs <- numpairs + 1
              if (letterij == letterkj) { numid <- numid + 1}
            }
          }
        }
        pcid <- (numid*100)/(numpairs)
        # Only consider this column if at least %minpcid of the pairs of letters are identical:
        if (pcid >= minpcid)
        {
          for (j in 1:numseqs)
          {
            seqj <- alignment$seq[[j]]
            letterij <- substr(seqj,i,i)
            newalignmentj <- newalignment$seq[[j]]
            newalignmentj <- paste(newalignmentj,letterij,sep="")
            newalignment$seq[[j]] <- newalignmentj
          }
        }
      }
    }
    return(newalignment)
  }
unrootedNJtree <- function(alignment,type, titlee){
  # this function requires the ape and seqinR packages:
  require("ape")
  require("seqinr")
  # define a function for making a tree:
  makemytree <- function(alignmentmat)
  {
    alignment <- ape::as.alignment(alignmentmat)
    if      (type == "protein")
    {
      mydist <- dist.alignment(alignment)
    }
    else if (type == "DNA")
    {
      alignmentbin <- as.DNAbin(alignment)
      mydist <- dist.dna(alignmentbin)
    }
    mytree <- nj(mydist)
    mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
    return(mytree)
  }
  # infer a tree
  mymat  <- as.matrix.alignment(alignment)
  mytree <- makemytree(mymat)
  # bootstrap the tree
  myboot <- boot.phylo(mytree, mymat, makemytree)
  # plot the tree:
  plot.phylo(mytree,type="u")   # plot the unrooted phylogenetic tree
  title(titlee)
  nodelabels(myboot,cex=0.7)    # plot the bootstrap values
  mytree$node.label <- myboot   # make the bootstrap values be the node labels
  return(mytree)
}
rootedNJtree <- function(alignment, theoutgroup, type, titlee){
  # load the ape and seqinR packages:
  require("ape")
  require("seqinr")
  # define a function for making a tree:
  makemytree <- function(alignmentmat, outgroup=`theoutgroup`, titl=titlee){
    # print(alignmentmat)
    alignment <- ape::as.alignment(alignmentmat)
    if (type == "protein")
    {
      mydist <- dist.alignment(alignment)
    }
    else if (type == "DNA")
    {
      alignmentbin <- as.DNAbin(alignment)
      mydist <- dist.dna(alignmentbin)
    }
    mytree <- nj(mydist)
    mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
    myrootedtree <- root(mytree, outgroup, r=TRUE)
    return(myrootedtree)
  }
  # infer a tree
  mymat  <- as.matrix.alignment(alignment)
  myrootedtree <- makemytree(mymat, outgroup=theoutgroup)
  # bootstrap the tree
  myboot <- boot.phylo(myrootedtree, mymat, makemytree)
  # plot the tree:
  plot.phylo(myrootedtree, type="p")  # plot the rooted phylogenetic tree
  title(titlee)
  nodelabels(myboot,cex=0.7)       # plot the bootstrap values
  myrootedtree$node.label <- myboot # make the bootstrap values be the node labels
  return(myrootedtree)
}

clean_alignment <- cleanAlignment(dirty_alignment, args$minpcnongap, args$minpcid)

n_str = grep(args$outgroup, clean_alignment$nam)
clean_alignment$nam[n_str]

if (args$type=='unrooted'){
  tree <- unrootedNJtree(clean_alignment, type=args$molecula)
  png(filename=paste("out", args$input, ".aln"))
  plot(tree)
  dev.off()
}

if (args$type=='rooted'){
  tree <- rootedNJtree(alignment=clean_alignment, theoutgroup=clean_alignment$nam[n_str], type=args$molecula, titlee=args$title)
  png(filename=paste("out/", args$input, ".png", sep = ''))
  plot(tree)
  dev.off()
}


  
