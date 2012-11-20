#! /usr/bin/python

import sys
# from progressbar import *
from BitVector import *

cladesAsMultifurcations = True

USAGE = "./script.py <ncbi.dmp> <fileWithRelevantTaxIds>"

if len(sys.argv) != 3:
    print (USAGE)
    sys.exit()

nodeDmp = sys.argv[1]
relevantTaxaFile = sys.argv[2]

def parseAdjList(nodeDmp):
    fh = open(nodeDmp)
    lines=fh.readlines()
    fh.close()

    adjList = {}
    for line in lines:
        tok = line.split('|')
        child = int(tok[0].strip())
        father = int(tok[1].strip())

        if(father == child):
            sys.stderr.write("excluding self-referential relationship of " + str(child) +  "\n" )
            continue
    
        if father in adjList:
            adjList[father].append(child)
        else:
            adjList[father] = [child]
    return adjList


def parseTaxList(relevantTaxaFile):
    fh = open(relevantTaxaFile)
    lines = fh.readlines()
    fh.close()
    return map(lambda x : int(x.strip()), lines)



def printSubtreeRecursive(adjList,taxaBv, currentTaxon):
    stringsOfChildren = []

    global found
    global i
    global numBits
    global unrooted
    
    i += 1
#     progress.update(i)

    if found < numBits and currentTaxon in adjList.keys():
        for child in adjList[currentTaxon]:
            string = printSubtreeRecursive(adjList, taxaBv, child)
            if string != "":
                stringsOfChildren.append(string)

    if cladesAsMultifurcations and taxaBv[currentTaxon] == 1 :
        found += 1
        stringsOfChildren.append(str(currentTaxon))

    result = ",".join(stringsOfChildren)

    if not cladesAsMultifurcations and taxaBv[currentTaxon] == 1:
        found += 1
        if result == "":
            result = str(currentTaxon)
        else:
            result = "(" +  result + ")" + str(currentTaxon)
    elif(len(stringsOfChildren) > 1):
        if(not unrooted):
            unrooted = True
        else:
            result = "(" + result + ")"

    return result


def printSubtree(adjList, taxaBv):
    result = printSubtreeRecursive(adjList, taxaBv, 1)
    sys.stdout.write( result + ";\n")

adjList = parseAdjList(nodeDmp)
relevantTaxa = parseTaxList(relevantTaxaFile)

maxi = 0
for elem in adjList.keys():
    tmp = max(adjList[elem])
    if tmp > maxi:
        maxi = tmp

taxaBv = BitVector(size=maxi)
for taxon in relevantTaxa:
    taxaBv[taxon] = 1
    
i = 0
found = 0
unrooted = False
numBits = taxaBv.count_bits()
# progress = ProgressBar(554623)

printSubtree(adjList, taxaBv)
