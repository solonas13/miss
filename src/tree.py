#! /usr/bin/python

from Bio import Phylo
from cStringIO import StringIO
import subprocess
import glob 
import random 
import os 
import re 

# This as a python class that contains both the tree with branch
# lengths and also computes the insertion lws for a given taxo. 
# 
# see the main function for an example




pattern = re.compile('^\s+"tree": "(.*)",\s*')
placePattern = re.compile('\s*{"p":\[(.*)\], "n":\[.*\]}')



# TODO 
raxml="/lhome/labererae/proj/mislabel-repo/lib/standard-RAxML/raxmlHPC-PTHREADS-SSE3"
threads="4"



def  maptotupleML(x) : 
    tok = x.split(",")
    return (int(tok[0]), float(tok[2]))


class Tree(object):
    def __init__(self, tree):
        self.treeFile = tree
        self.origTree = Phylo.read(tree, "newick")
        self.lwTree = None

    def computeLWs(self, alnFile, ml=True): 
        method = "v"
        if not ml: 
            method = "y"
            
        theid = str(random.randrange(999999999))
        run = [raxml , 
               "-T", threads,
               "-s", alnFile,
               "-t", self.treeFile,
               "-m",  "GTRGAMMA",
               "-f", method,
               "-n" , theid 
               ]
        
        with open(os.devnull, "w") as fnull: 
            bla = subprocess.call(run, stdout = fnull)
        jsonFile = glob.glob(os.getcwd() + "/RAxML*." + theid + ".jplace")[0]        
        string = self.__createTreeFileWithLws(jsonFile, ml)
        self.lwTree = Phylo.read(StringIO(string), "newick") 
        c = self.lwTree.find_clades().next()
        c.branch_length=0
        
        # cleanup 
        for elem in glob.glob(os.getcwd() + "/*." + theid): 
            os.remove(elem)    
        for elem in glob.glob(os.getcwd() + "/RAxML*." + theid + ".jplace"): 
            os.remove(elem)



    def __createTreeFileWithLws(self, jsonFile, ml=True ): 
        lines = open(jsonFile, "r").readlines()

        treeMatches = filter(lambda x : x.strip().startswith('"tree":') ,lines)
        assert(len(treeMatches) == 1)
        treeMatches = treeMatches[0]
        match = pattern.match(treeMatches)
        assert(match)        
        tree = match.group(1)

        result = ""
        parseNow = False
        for line in lines: 
            if parseNow : 
                match = placePattern.match(line.strip())
                assert(match)
                result = match.group(1)
                break
            else: 
                if(line.strip().startswith('"placements"')): 
                    parseNow = True
        
        assert(result)

        placements = filter(lambda x : x != "" , map(lambda x :  x.strip("[,") , result.split("]")))

        if ml: 
            placeDict = {}
            for p in map( maptotupleML, placements) : 
                placeDict[p[0]] = p[1]        
            tree = re.sub(":[0-9\.]*{([0-9]*)}", 
                          lambda x : ':' + str(placeDict[int(x.group(1))]),tree ) 
        else : 
            placeDict = {} 
            num = float(len(placements)) 
            for p in placements: 
                tok = p.split(",")
                placeDict[int(tok[0])] = float(1) / num
            tree = re.sub(":[0-9\.]*{([0-9]*)}", 
                          lambda x :   ":" + str(placeDict.get(int(x.group(1)), 0.0)) , tree ) 
        return  tree


    def __str__(self):
        return str(self.origTree) + "\n" + str(self.lwTree)

    
    def scoreBipartition(self, bip, withComplement = False): 
        "Uses the likelihood weights in to score each bipartition in the tree using bip as reference."
        allLeaves = set(map(lambda x : x.name ,  self.lwTree.get_terminals() ) )
        bip = set(bip)

        score = 0
        
        for clade in self.lwTree.find_clades(): 
            lw = clade.branch_length
            
            bipHere = set(map(lambda x : x.name , clade.get_terminals()))             
            complBipHere = allLeaves - bipHere 

            scoreA = float(len(bip & bipHere) ) / float(len(bipHere)) 
            if withComplement and len(complBipHere) != 0 : 
                scoreB = float(len(bip & complBipHere) ) / float(len(complBipHere)) 
                scoreA = max(scoreA,scoreB)
            score += scoreA * lw 
        print score 



if __name__ == "__main__": 
    treeFile = "/lhome/labererae/proj/mislabel-repo/data/trees/150.tre"
    alnFile = "/lhome/labererae/proj/mislabel-repo/data/aln/150.aln"
    prunedTree = "prunedTree.tre" 

    tree = Phylo.read(treeFile, "newick")
    tree.prune("Species144")
    Phylo.write(tree, prunedTree, "newick")

    lwTree = Tree(prunedTree)
    lwTree.computeLWs(alnFile, True)
    os.remove(prunedTree)

    print lwTree.lwTree.format("newick")
    print lwTree.scoreBipartition(["Species005", "Species202", "Species026"])
