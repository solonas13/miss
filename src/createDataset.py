#! /usr/bin/python

from tree import *
import random 
from Bio import Phylo
from cStringIO import StringIO


USAGE ="""
./script <aln> <tree>

"""


import sys

if len(sys.argv) != 3: 
   print USAGE
   sys.exit()


alnFile=sys.argv[1]
treeFile=sys.argv[2]

tree = Phylo.read(treeFile, "newick")
taxa = map(lambda x : x.name ,  tree.get_terminals())
toMislabel = random.choice(taxa) 
print "mislabelling "  + toMislabel


## :TODO: we must ensure in tree.py that the tree is unrooted ... => use newick utils 

tok = sys.argv[0].split("/") 
tok = tok[0:-1]
rootDir =  "/".join(tok)
reroot_cmd=rootDir + "../lib/newick-utils-1.6/src/nw_reroot"



lwTreeByTaxon = {}
for taxon in taxa  : 
    print "evaluating lw for " + taxon
    tmpFile = "pruned.tre"
    tree = Phylo.read(treeFile, "newick")
    tree.prune(taxon)
    Phylo.write(tree, tmpFile, "newick")
    
    lwTree = Tree(tmpFile)
    lwTree.computeLWs(alnFile, False)
    lwTreeByTaxon[taxon] = lwTree
    os.remove(tmpFile)


