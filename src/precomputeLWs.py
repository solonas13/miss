#! /usr/bin/python


# This is a very small script, that can be used to compute the
# likelihood weights on the server. For the evaluation step, this is
# important, we do not want to recompute lws each time with ML this is
# terribly expensive. 



import sys 
from tree import * 
from taxonomy import * 


if len(sys.argv) != 5 : 
    print "USAGE: ./script <treeFile> <alnFile> <id> <useML>"
    sys.exit()

treeFile = sys.argv[1]
alnFile = sys.argv[2] 
theid = sys.argv[3]
useML = sys.argv[4] == "True"

tree = LwTree(treeFile)

leaves = tree.get_leaf_names()

for taxon in leaves: 
    tr = LwTree(treeFile)
    tr.pruneTaxa([taxon])
    tr.computeLWs(alnFile, useML)
    tr.write(features=["lw"], format=2, outfile=theid + "." + taxon + "."  + ("pars" if not useML else "ml" )+ ".tre")
    
