#! /usr/bin/python

import sys 
import random 
from copy import * 
from taxonomy import * 


if len(sys.argv) != 7 : 
    print "./script <taxonomy> <givenTaxa> <renamingFile> <numOfMisLabel> <depthOfMislabeling> <runid>"
    sys.exit()

numMislabel = int(sys.argv[4])
depthOfMislabeling = int(sys.argv[5])
runId = sys.argv[6]

taxonomyFn = sys.argv[1]

tax = Taxonomy()
tax.init_parseTaxFile(taxonomyFn)

relevantTaxa  = map(lambda x: x.strip(), open(sys.argv[2], "r").readlines()) 


renamings = map(lambda x : x.strip(), open(sys.argv[3], "r").readlines())
renamDict = {}
for elem in renamings: 
    [origValue, newName] = elem.split()
    renamDict[origValue]  = newName

reverseDict = {}
for elem in renamings:
    [origValue,newName] = elem.split()
    reverseDict[newName] =origValue

renamedRelevant = []
for elem in relevantTaxa: 
    if renamDict.has_key(elem):
        renamedRelevant.append(renamDict[elem])
    else  : 
        renamedRelevant.append(elem)

tax.reduceTaxonomy(renamedRelevant)
correctTaxanomy = deepcopy(tax)

leaves = list(tax.getLeaves())
random.shuffle(leaves)

ctr = 0
mislabeles = []
correctTaxonomy = deepcopy(tax)
oldTax = deepcopy(tax)
mislabeled = 0
while(mislabeled < numMislabel): 
    currentTax = leaves[ctr]
    tax.mislabel(currentTax, depthOfMislabeling,depthOfMislabeling) 
    
    if not (oldTax == tax)  : 
        mislabeles.append(currentTax)
        mislabeled += 1 
        oldTax = deepcopy(tax)
    else : 
        tax = deepcopy(oldTax)
    ctr += 1 


correctTaxonomy.renameLeaves(reverseDict)
correctTaxonomy.saveToFile("correctTaxonomy."  + runId)
tax.renameLeaves(reverseDict)
tax.saveToFile("mislabeledTaxonomy." + runId)

fh = open("info." + runId, "w")
fh.write("mislabeled=%s\n" % ",".join(mislabeles))
fh.write("attemptedLevel=%s\n" % depthOfMislabeling)
fh.write("correctNewick=%s\n" % correctTaxonomy.getNewickString() ) 
fh.write("mislabeledNewick=%s\n"  % tax.getNewickString() ) 
fh.close()


