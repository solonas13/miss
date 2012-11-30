#! /usr/bin/python

# So this an example file of how to score specific taxa on
# mislabelledness.

import sys 
from tree import * 
from taxonomy import * 


USAGE= """./script <treeFile> <alnFile> <taxonomy> <useML> <infoFile>
* treefile: a ML tree 
* aln file: the alignment file 
* taxonomy: a taxonomy file (=> leaves in the tax must be the same as in the tree!)
* useML: True, if ML should be used for placements, False for parsimony
* infoFile: the info file, specifying which taxon has been mislabelled.  
"""



if __name__ == '__main__': 

    if len(sys.argv)  != 6: 
        print USAGE
        sys.exit()

    basePath = os.path.dirname(sys.argv[0]) 
    treeFile = sys.argv[1]
    alnFile = sys.argv[2]
    taxFile = sys.argv[3]
    useML = sys.argv[4] == "True"
    infoFile = sys.argv[5]
    
    # parse the taxonomy 
    tax = Taxonomy()
    tax.init_parseTaxFile(taxFile)

    # get the mislabelled taxa 
    mislabeled = set(open(infoFile, "r").readline().strip().split("=")[1].split(",")) 
    
    # evaluate the taxa 
    # print "taxon\tlcaScore\toverlapScore"
    leaves = tax.getLeaves()
    for taxon in leaves:         
        potentialPrecomputedFile = os.path.dirname(alnFile) + "/../precomputedLWs/" + os.path.splitext(os.path.basename(alnFile))[0] + "." + taxon + '.' + ( "ml" if useML else "pars" ) + ".tre"
        if(os.path.exists(potentialPrecomputedFile)): 
            sys.stderr.write("using file %s\n" % potentialPrecomputedFile)
            tr = LwTree(potentialPrecomputedFile) 
            tr.regenerateLws()
        else :             
            sys.stderr.write("not using file %s\n" % potentialPrecomputedFile) 
            tr = LwTree(treeFile)
            tr.pruneTaxa([taxon])
            tr.computeLWs(alnFile, useML)

        # this is stupid =/ 
        bip = set()
        for elem in tax.getNthBipartition(taxon,1): 
            if elem != taxon : 
                bip.add(elem)
        assert(not (taxon in bip)) 

        lcaScore = tr.getLcaScore(bip)
        overlapScore = tr.getOverlapScore(bip)
        
        print "%s\t%f\t%f\t%s" % (taxon,  lcaScore, overlapScore, ( "mislabel" if (taxon in mislabeled)  else "correct" ) )

        
