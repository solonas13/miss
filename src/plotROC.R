#! /usr/bin/Rscript

# This script takes the output of various runs of evaluateADataset.py
# and plots a ROC.

library(verification)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)  != 2 )
  {
    cat("Usage: ./script <id> <file>\nwhere <file> is the output of evaluate.py (best: multiple runs). \n<id> serves as a name for the plot.\n")
    q()
  }

tab = read.table(args[2])

reality = as.character(tab[,4]) 
reality[reality == "mislabel"] = 0
reality[reality == "correct"] = 1
reality = as.numeric(reality) 
predictions = tab[,c(2,3)]

pdf(paste(args[1], ".pdf", sep=""))
roc.plot(reality,predictions, legend=T, leg.text = c("LCA", "overlapScore"))
bla = dev.off()
