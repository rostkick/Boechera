#!/bin/bash
# merge eggnog's proteins with won proteins
# for merging you need patterns.txt file
paste <(find . -name *.fasta | sort -k2 -t"/") <(find . -name *.fa | grep -f patterns.txt | sort -k2 -t"/" | grep -v FIS) | while read i j; do cat $i $j > out/${j##*/}; done

