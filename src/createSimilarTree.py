#! /usr/bin/python

from tempfile import * 
from random import *
from ete2 import * 
import sys 

raxml="/lhome/labererae/proj/standard-RAxML/raxmlHPC-SSE3"
raxmlPara="/lhome/labererae/proj/standard-RAxML/raxmlHPC-PTHREADS-SSE3"
prs ="/lhome/labererae/proj/Parsimonator-1.0.2/parsimonator-SSE3"
numCores = 3 

def computeParsimonyTrees(treeFile,alnFile, runId, num=100):
    """ Computes a bunch of parsimony trees and returns the base file name  
    """

    seed = randrange(0,99999999999)
    prscall = " ".join([
            prs, 
            "-n", runId, 
            "-N", str(num),
            "-s" , alnFile,
            "-t", treeFile,
            "-p", str(seed) ,
            " > /dev/null "
            ])
    os.system(prscall) 
    os.system("rm %s" % "RAxML_info." + runId )
    return "RAxML_parsimonyTree." + runId + "."
    


def getOptimizedTree(treeFile, alnFile, runId):
    nniCall = " ".join([            
            raxmlPara, 
            "-T", str(numCores),
            "-n nni." + runId,
            "-m", "GTRGAMMA", 
            "-f", "J", 
            "-s", alnFile,
            "-t", treeFile
            ])
    os.system(nniCall + " > /dev/null ")
    result = "RAxML_fastTreeSH_Support.nni."  + runId

    fh = open("RAxML_info.nni." + runId , "r")
    lnl = filter(lambda x : x.startswith("Final Likelihood of NNI-optimized"), fh.readlines())
    fh.close() 
    assert(len(lnl) == 1 )
    likelihood = float(lnl[0].split(":")[1].strip())
    
    fh = open(result )
    tree = fh.readlines()
    fh.close()
    assert(len(tree) == 1 )
    tree = tree[0].strip()

    os.system("rm %s " % treeFile)
    # os.system("rm RAxML_*.eval.%s"  % runId)
    os.system("rm RAxML_*.nni.%s"  % runId)

    return [ tree, likelihood ] 




def getrelRFDistance(origTreeFileName, tree, runId): 
    tmpFileName = mkstemp()[1]
    
    fh = open(origTreeFileName, "r")
    origTree = fh.readline().strip()
    fh.close()   
    
    fh = open(tmpFileName, "w")
    fh.write(origTree + "\n")
    fh.write(tree + "\n" )
    fh.close()

    call = " ".join([
            raxml,             
            "-m", "GTRCAT", 
            "-n", "rf." + runId ,
            "-f", "r",
            "-z", tmpFileName
            ]) 
    os.system(call + " > /dev/null")

    os.system("rm %s " % tmpFileName)
    fh = open("RAxML_RF-Distances.rf." +runId, "r")
    relRF = float(fh.readline().split(" ")[3])
    fh.close()
    os.system("rm RAxML_*." + "rf." + runId)

    return relRF




def getSimilarTrees(alnFile,treeFile, runId, proportion = 0.25, num =100 , minRF=0.03, maxRF=0.1):
    tr = Tree(treeFile)
    origLeaves = tr.get_leaf_names()
    alist = origLeaves
    shuffle(alist)
    toPrune = alist[0:int(len(origLeaves) * proportion)]

    toStay = set(origLeaves) -  set(toPrune)
    tr.prune(toStay)

    assert(len(tr) + len(toPrune) == len(origLeaves))
    tr.unroot()
    
    
    tmp = mkstemp()
    tr.write(format=2, outfile=tmp[1])
    
    baseName = computeParsimonyTrees(tmp[1], alnFile, runId, num=num)
    resultList = []
    for i in range(0,num): 
        result = getOptimizedTree("RAxML_parsimonyTree.%s.%d" % ( runId , i), alnFile, runId)
        result.append(getrelRFDistance(treeFile, result[0], runId))
        resultList.append(result) 
        sys.stderr.write( "%f\t%f\n" % (result[1], result[2]))
    os.system("rm %s"  % tmp[1])

    result = filter(lambda x : minRF <= x[2] and x[2] <= maxRF, resultList )
    tooSmall = len(filter(lambda x : x[2]  < minRF ,  resultList))
    tooLarge = len(filter(lambda x : maxRF < x[2]  ,  resultList))
    return (result, tooSmall, tooLarge )



def floatEql(a,b): 
    return  abs(a - b)  < 0.00001

def setify(alist): 
    result = [] 
    for elem in alist: 
        isNew = True 
        for alreadyThere in result: 
            equal = floatEql(elem[1], alreadyThere[1]) and floatEql(elem[2], alreadyThere[2])
            if equal : 
                isNew =  False 
        if isNew: 
            result.append(elem)
    return result




if len(sys.argv ) != 7 : 
    print """./script <aln> <tree> <num> <rfMin> <rfMax> <id>
* aln: alignment file 
* tree: tree file from which trees will be derived
* num: number of output trees 
* rfMin: minimum rf-distance
* rfMax: maximum rf-distance
* id : an id for this run
"""
    sys.exit()


alnFile = sys.argv[1]
treeFile = sys.argv[2]
num = int(sys.argv[3])
minRF= float(sys.argv[4])
maxRF= float(sys.argv[5])
theId = sys.argv[6]

assert(minRF < maxRF)

trees = []
currentProportion = 0.2
step = 0.025
while(len(trees) < num): 
    sys.stderr.write("numTrees = %d/%d\t current proportion= %f\n" % (len(trees), num,   currentProportion))
    remaining = num - len(trees)
    tryingThisRound = min(5,remaining)
    (treesHere, tooSmall, tooLarge) = getSimilarTrees(alnFile, treeFile, theId, proportion=currentProportion, num= tryingThisRound, minRF=minRF, maxRF=maxRF) 
    trees.extend(treesHere)

    # modify the cutoff 
    if(tryingThisRound >= 3 and tooLarge + tooSmall > tryingThisRound / 2 ): 
        if tooSmall > tooLarge : 
            currentProportion += step
        else : 
            currentProportion -= step
        if currentProportion <= 0: 
            currentProportion = step

    # remove identical 
    trees = setify(trees)

fh = open("simTree.result." + theId, "w")
infoFH = open("simTree.info." + theId, "w") 
for elem in trees: 
    fh.write("\n".join(map(lambda x : x[0], trees)) + "\n")
    infoFH.write("\n".join(map(lambda x : "\t".join([str(x[1]),str(x[2])]), trees)) + "\n")

fh.close()
infoFH.close()
