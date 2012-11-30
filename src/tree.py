#! /usr/bin/python 


import os
import sys 

import json
import glob 
import re 
import random
import subprocess

from ete2 import * 
from pprint import pprint
from tempfile import  * 


numCores = 4 
raxmlPath = os.path.dirname(sys.argv[0]) + "/../lib/standard-RAxML/raxmlHPC-PTHREADS-SSE3"


class LwTree(Tree):     
    """ A tree derived from the ete2 tree including likelihood weights """ 

    
    # TODO does not work =( 
    # def __init__(self, treeFile): 
    #     if not os.path.exists(raxmlPath): 
    #         sys.stderr.write("expecting path to raxml: %s\nPlease install!" % raxmlPath)
    #     Tree.init_self
        # super(Tree, self).__init__(treeFile)
        # if hasattr(self.get_descendants()[0], 'lw'): 
        #     print "has it"
        #     for desc in self.get_descendants():
        #         desc.lw = float(desc.lw)
        # else : 
        #     print "has not it "
    
    def regenerateLws(self): 
        """
        This is a hack. If we re-read annotated trees (because we
        precomputed the likelihood weights), then the labels are
        incorrectly parsed as strings instead of floats.
        """
        for desc in self.get_descendants():
            desc.lw = float(desc.lw)

                
            
            
    def pruneTaxa(self, listOfTaxa):         
        "Prunes a list of taxa. This inverts the method provided by ete2 for convenience. "
        self.prune(list(set(self.get_leaf_names()) - set(listOfTaxa)))



    def __annotate__(self, lwTree): 
        assert( 1.0 - sum(map(lambda x : x.dist, lwTree.get_descendants()))  < 0.0001 )

        for desc in lwTree.get_descendants(): 
            if desc.is_leaf(): 
                res = self.search_nodes(name=desc.get_leaf_names()[0])
                assert(len(res) == 1 )
                res[0].add_feature("lw", desc.dist) 
            else : 
                leaves = desc.get_leaf_names()
                ancNode = self.get_common_ancestor(leaves)
                assert(len(leaves) == len(ancNode.get_leaf_names()))
                ancNode.add_feature("lw",  desc.dist) 

        # did we get all lws? 
        assert(1.0 - sum(  map(lambda x : x.lw, self.get_descendants()))  < 0.0001  ) 


    def __annotateWithParsLws(self, jsonString): 
        data =  json.loads(jsonString)

        positions = map(lambda x : x[0], data['placements'][0]['p']) 
        lws = 1 / float(len(positions)) 
        
        lwMap = {}
        for i in range(0, len(positions)): 
            lwMap[positions[i]] = lws

        treeString = data['tree'] 
        treeString = re.sub(":[0-9\.]*{([0-9]*)}",
                            lambda x : 
                            ((":" +  str(lwMap[int(x.group(1))]) ) 
                             if  lwMap.has_key(int(x.group(1))) 
                             else ":0.0" ) 
                            ,treeString ) 
        lwTree = Tree(treeString, format=1)
        self.__annotate__(lwTree)


    def __annotateWithAllLws(self, jsonString):
        data =  json.loads(jsonString)

        positions = map(lambda x : x[0], data['placements'][0]['p']) 
        lws = map(lambda x : x[2], data['placements'][0]['p']) 

        lwMap = {}
        for i in range(0, len(positions)): 
            lwMap[positions[i]] = lws[i]

        treeString = data['tree']
        
        treeString = re.sub(":[0-9\.]*{([0-9]*)}", 
                      lambda x :
                          ':' + str(lwMap[int(x.group(1))]),treeString ) 
        lwTree = Tree(treeString, format=1)
        self.__annotate__(lwTree)



    def computeLWs(self, alnFile, ml=True): 
        """ Computes likelihood weights and annotates it to the
        tree. Important: you'll have to assure yourself, that the
        alignment only contains one more taxon than the tree (e.g.,
        prune first)"""         
        
        tmpFile = mkstemp()
        tmpFileName = tmpFile[1]
        
        self.unroot()
        self.write(format=1, outfile=tmpFileName)
        
        # prepare raxml run 
        theid = str(random.randrange(999999999))

        run = [raxmlPath , 
               "-T", str(numCores),
               "-s", alnFile,
               "-t", tmpFileName,
               "-m",  "GTRGAMMA",
               "-f", ("v" if ml else "y" ),
               "-n" , theid 
               ]

        # :TRICKY: if a partition file is needed, we'll determine that
        # here. Notice, that this requires, that the partition file
        # bears the same name the aln file, but with file extension
        # .par.         
        potentialPartitionFile = os.path.splitext(alnFile)[0] + ".par"        
        if os.path.exists(potentialPartitionFile): 
            run.append("-q")
            run.append(potentialPartitionFile)
        
        with open(os.devnull, "w") as fnull: 
            bla = subprocess.call(run, stdout = fnull)
        os.remove(tmpFileName)

        if(ml): 
            self.__annotateWithAllLws(open("RAxML_portableTree." + theid + ".jplace", "r").read())
        else : 
            self.__annotateWithParsLws(open("RAxML_portableTree." + theid + ".jplace", "r").read())

        # cleanup 
        for elem in glob.glob(os.getcwd() + "/*." + theid): 
            os.remove(elem)    
        for elem in glob.glob(os.getcwd() + "/RAxML*." + theid + ".jplace"): 
            os.remove(elem)


    def getLcaScore(self, bipartition): 
        """
        The lca score: find the lca of all the taxa in bipartition and
        sum up the liklihood weights in the subtree below the lca. 
        Thus the score ranges in [0,1]
        """ 
        # :TODO: unrooting necessary? 
        # self.unroot()

        bipartition = set(bipartition)

        # if len(bipartition) == 1 :
        #     res = self.search_nodes(name=desc.get_leaf_names()[0])
        #     assert(len(res) == 1) 
        #     score =  res.lw 
        # else : 
        subtree = self.get_common_ancestor(list(bipartition))
        score = sum(map(lambda x : x.lw, subtree.get_descendants())) 

        return score 


    def getOverlapScore(self, bipartition): 
        """
        The score, Andre proposed: the likelihood weight at each
        branch is weighted by the proportion of leaves that are of the
        same rank as the mislabel candidate.
        """
        # :TODO: unrooting necessary? 
        # self.unroot()
        
        bipartition = set(list(bipartition)) 
        # bipComplement = set(self.get_leaf_names()) - bipartition

        score =  sum(map(lambda x :  x.lw *  float(len(bipartition &  set(x.get_leaf_names() ) ))  /  float(len(set(x.get_leaf_names() ) )) ,
                         self.get_descendants()) ) 
        
        # sometimes summing leads to values > 1.0
        return min(score,1.0) 
        


USAGE= """./script <tree> <taxon> <aln> 

Compute likelihood weights for <taxon> on <tree> (still including
<taxon>). The reference alignment <aln> also already includes the
taxon.
"""

if __name__ == "__main__": 
    if len(sys.argv) != 4 : 
        print USAGE
        sys.exit()

    tree = sys.argv[1]
    taxon = sys.argv[2]
    aln = sys.argv[3]

    t = LwTree(tree)
    t.pruneTaxa([taxon])
    t.computeLWs(aln, False)
    
    print t.write(format=2, features = ["lw"])


    
