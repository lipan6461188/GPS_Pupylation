#!/bin/bash

perl covertSeq.pl --input 268Sequences --output sequence

perl fetchSubSeq.pl --upper 7 --down 7 --patter K --delimit 1 --input "sequence" --output "samples" --positive "POSITIVE.txt"

perl gps.pl --ubound 4 --dbound 0 --gap 0.001 --blosum BLOSUM62 --sample samples --output ROC --silent

Rscript plot.R