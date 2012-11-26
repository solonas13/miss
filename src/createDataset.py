#! /usr/bin/python




USAGE ="""./script <id> <aln> <numMisslabel> <mislabelDepth> [<tree>] 

* aln: alignment file 
* numMisslabel: number of mislabels to create in the dataset 
* mislabelDepth: number of ranks by which the taxon gets mislabelled. 
* tree: ml tree. If option is provided, a taxonomy will be created
  from the tree, otherwise a parsimony tree will be created from the
  alignment.
"""

numThreads = 4
maxiter = 10 
maxDist = 0.1 


import sys
import os 
import random 
from taxonomy import * 

import copy


############################
# parse command line args  #
############################
raxml=os.path.dirname( sys.argv[0])  + "/../lib/standard-RAxML/raxmlHPC-PTHREADS-SSE3"

try:
   with open(raxml) as f: pass
except IOError as e:
   print "Please install standard raxml into the lib folder and compile the non-pthreads sse3-version."
   sys.exit()


if len(sys.argv) < 4 : 
   print USAGE
   sys.exit()

theId = sys.argv[1]
alnFile=sys.argv[2]
numMiss = int(sys.argv[3]) 
mislabelRank = int(sys.argv[4])

createParsTree = True
treeFile = ""
if len(sys.argv) == 6: 
   createParsTree = False
   treeFile = sys.argv[5]



##############
# functions  #
##############
def createSimilarTree(alnFile, maxDist, maxIter, theId, useParsimony=True): 
   """ 
   create a parsimony tree from the alignment file <alnFile> with
   relative RF-distance < <maxDist> (float).
   The raxml id for this run is <theId>. 

   If you specify useParsimony as False, then instead of a parsimony,
   we will create a quick&dirty fastTree with raxml -f F. 
   
   The program will give up after <maxIter> iterations. Then, it
   yields the closest tree encountered so far.

   If a partition file exists in the same folder, then it will be used
   (important for protein data).
   """ 
   distanceTooHigh = True
   seed = random.randint(0,1000000000)

   baseName = os.path.splitext(alnFile)[0]
   if os.path.exists(baseName + ".par") : 
      partitionString = "-q " +  baseName + ".par"
   else : 
      partitionString = ""


   origTree = baseName + ".tre"
   if not os.path.exists(origTree): 
      print "Expected tree " + origTree + 'to exist. It does not.' 
      sys.exit()

   if useParsimony: 
      resultBase = "RAxML_parsimonyTree."
   else : 
      resultBase = "RAxML_fastTree." 
   resultTree = resultBase + theId

   bestDistance = 1.0
   bestTree = resultBase + theId + ".saved"
   while distanceTooHigh: 
      parsimonyString = "-y" if useParsimony else "-f E"

      treeCmd = " ".join([
            raxml, 
            "-T", str(numThreads), 
            "-s" , alnFile, 
            partitionString, 
            "-n" , theId, 
            parsimonyString, 
            "-m", "GTRCAT",
            "-p" , str(seed),
            "> /dev/null"
            ])  

      os.system(treeCmd)

      treeFile = resultBase + theId
      cmd = "cat " + treeFile + " " + origTree + " > bothTrees"
      os.system(cmd) 

      # compute RF distance 
      rfCmd = " ".join([
            raxml, 
            "-T", str(numThreads), 
            "-s" , alnFile, 
            partitionString, 
            "-n" , theId + "Eval", 
            "-m", "GTRCAT",
            "-z", "bothTrees",
            "-f r > /dev/null", 
            ])
      os.system(rfCmd)
      
      rfFile = "RAxML_RF-Distances." + theId + "Eval"
      rfDistance = float(open(rfFile, "r").readline().strip().split(" ")[3])
      print "The rf-distance was " + str(rfDistance)

      if rfDistance < bestDistance and rfDistance != 0: 
         os.system("cp " + treeFile + " "  +  bestTree) 
         bestDistance = rfDistance 

      os.system("rm bothTrees *." + theId + " *." + theId + "Eval")
      seed += 1 
      distanceTooHigh =  rfDistance > maxDist or rfDistance == 0

   return (bestDistance, bestTree)

#########
# main  #
#########
if __name__ == "__main__": 
   refTree = ""
   if createParsTree:       
      result = createSimilarTree(alnFile, maxDist, maxiter, str(random.randint(0,99999999999)), False)
      refTree = result[1]
   else : 
      refTree = treeFile
      result = (0, "")
   
   tax = Taxonomy()
   tax.init_extractRandomlyFromTree(refTree)

   initTax = copy.deepcopy(tax)

   # mislabeling the taxa 
   leaves = list(tax.getLeaves())
   random.shuffle(leaves)   
   toMislabel = leaves
   mislabeled = []

   for toMis in  toMislabel: 
      currentTax = copy.deepcopy(tax)
      ctr = 0
      while tax == currentTax: 
         tax.mislabel(toMis, mislabelRank, mislabelRank)
         if ctr == 5 : 
            break
         ctr += 1 
 
      if not ( tax ==  currentTax )  : 
         mislabeled.append(toMis)
      if len(mislabeled) == numMiss: 
         break 

   # print the results 
   initTax.saveToFile("correctTaxonomy." + theId + ".tax") 
   tax.saveToFile("mislabeledTaxonomy." + theId + ".tax")

   # write info file 
   fh = open("info." + theId + ".txt", "w")
   fh.write("mislabelledTaxa=%s\n" % ",".join(mislabeled))   
   fh.write("distanceOfTaxonomyToRefTree=%f\n" % result[0])
   fh.write("mislabelLevel=%d\n" % (mislabelRank-1))
   fh.write("initialAlignment=%s\n" % alnFile)
   fh.write("correctTaxonomy=%s\n" % initTax.getNewickString())
   fh.write("mislabelledTaxonomy=%s\n" % tax.getNewickString())
   fh.close()
