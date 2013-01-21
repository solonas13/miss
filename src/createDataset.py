#! /usr/bin/python

import sys
import random
from copy import * 
from taxonomy  import * 


if len(sys.argv) != 9 : 
    print """ USAGE: ./script <tree> <simTree>  <numLevels> <patDist> <midPoint> <numMisLabel> <ranksUp> <runId> 
with
* tree 
* numLevels: depth of the taxonomy 
* simTree: either a similar tree or or the same tree the topology was extracted from (for checking later)
* patDist: True, if patristic distance should be used (otherwise non-patristics)  
* midPoint: True, if taxonomy should be rooted by midpoint, otherwise innermost node is used for rooting
* numMisLabel: number of mislabelled taxa 
* ranksUp: indicate by how many ranks taxa are mislabelled
* runId: runid 
"""
    sys.exit()

treeFile =  sys.argv[1]
simTree = sys.argv[2]
numLevels = int(sys.argv[3])
patDist = sys.argv[4] ==  "True"
midPoint = sys.argv[5] == "True"
numMislabel = int(sys.argv[6]) 
ranksUp = int(sys.argv[7])
runId = sys.argv[8]


needInconsistency = treeFile != simTree

if not needInconsistency: 
    tax = Taxonomy()
    nameMapping = tax.init_extractRandomlyFromTreeImproved(simTree, maxlevel = numLevels, usePatristicDistance = patDist, useMidpoint = midPoint)
    assert(tax.taxonomyFitsTree(treeFile,nameMapping))
else : 
    treeIsConsistentWithTaxonom = True 
    while treeIsConsistentWithTaxonom: 
        tax = Taxonomy()
        nameMapping = tax.init_extractRandomlyFromTreeImproved(simTree, maxlevel = numLevels, usePatristicDistance = patDist, useMidpoint = midPoint)
        treeIsConsistentWithTaxonom = tax.taxonomyFitsTree(treeFile,nameMapping)
    
leaves = list(tax.getLeaves())
random.shuffle(leaves)

ctr = 0
mislabeles = []
correctTaxonomy = deepcopy(tax)
oldTax = deepcopy(tax)
mislabeled = 0
while(mislabeled < numMislabel): 
    currentTax = leaves[ctr]
    tax.mislabel(currentTax, ranksUp, ranksUp) 
    
    if not (oldTax == tax)  : 
        mislabeles.append(currentTax)
        mislabeled += 1 
        oldTax = deepcopy(tax)
    else : 
        tax = deepcopy(oldTax)
    ctr += 1 



# write dataset 
fh = open("nameMapping." + runId, "w")
fh.write( "\n".join(map(lambda x : x[0] + "\t" + x[1] ,  nameMapping.items())) + "\n")
fh.close()

correctTaxonomy.saveToFile("correctTaxonomy."  + runId)
tax.saveToFile("mislabeledTaxonomy." + runId)

fh = open("info." + runId, "w")
fh.write("mislabeled=%s\n" % ",".join(mislabeles))
fh.write("treeFile=%s\n" % treeFile)
fh.write("similarTree=%s\n" % simTree)
fh.write("patDist=%s\n" % str(patDist))
fh.write("midPoint=%s\n" % str(midPoint))
fh.write("ranksUp=%s\n" % ranksUp)
fh.close()

