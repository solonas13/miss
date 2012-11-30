#! /usr/bin/python

import sys 
import random 
from ete2 import *
import copy 


# a flexible representation of a taxonomy that
# * allows to reduce the taxonomy to some specified taxa (much faster than previously)
# * creates a random taxonomy from a tree 
# * mislabels a taxon in a taxonmy (i.e., the taxon is moved around in the taxonomy)

# => see the main function for an example

# known problems: internal nodes that are a species (and additionally
# have a sub-species) => but this should not occur unless you have a
# really crazy dataset.




class IdFactory(object): 
    "A small helper class that provides ids, when a taxonomy is created from a tree."
    def __init__(self, leaves, start): 
        self.blocked = set(leaves)
        self.nextId = start

    def produceId(self): 
        assert(str(self.nextId) not in self.blocked)
        result = self.nextId        
        self.nextId += 1
        while(str(self.nextId) in self.blocked): 
            self.nextId += 1 
        return str(result)
        
        


class Taxonomy(object): 
    """ 
    Represents a taxonomy. Create an object and then either parse an
    existing taxonomy file or use a raxml-style tree and reduce it
    randomly to a taxonomy that is consistent with this tree.

    For parsing the ncbi taxonomy, you have to translate the ncbi
    taxonomy into a 2-column format first. Something like 
    $ cat  nodes.dmp  | cut -f 1,2  -d '|' | tr -d  "\t"  | tr "|" "\\t" > convertedNcbi.tab
    should do. 
    
    """ 



    def __init__(self): 
        "Note: when parsing taxonomies with inner nodes, you may want to cleanup the taxonomy right after parsing to throw out inner nodes with only one child. "
        self.root = None
        self.childToParent = {}
        self.parentToChildren = {}
        self.dirty = True    


    def init_parseTaxFile(self, fn): 
        "parse an existing taxonomy from file <fn>" 
        lines = open(fn, "r").readlines()
        self.dirty = True
        for line in lines: 
            tok = line.strip().split()

            child = tok[0]
            parent = tok[1]

            if(child == parent): 
                assert(self.root == None)
                self.root = child 
                continue 
            
            assert(not self.childToParent.has_key(child))
            self.childToParent[child] = parent
            if self.parentToChildren.has_key(parent): 
                self.parentToChildren[parent].append(child)
            else : 
                self.parentToChildren[parent] = [child]

        if self.root == None: 
            self.root = self.__regenerateRoot()
            # print "determined root to be " + self.root
            
    def __regenerateRoot(self):
        leaves = self.getLeaves()
        innerNodes = set(self.parentToChildren.keys()) 
        innerNodesWoRoot = set(self.childToParent.keys())  - leaves
        assert(len(innerNodes) == len(innerNodesWoRoot) + 1)
        root = innerNodes - innerNodesWoRoot
        root = list(root)[0]
        return root


    def saveToFile(self, fn): 
        "Saves the taxonomy to a file <fn>. This file can later be read again with the function init_parseTaxFile()."
        fh = open(fn, "w")
        for child in self.childToParent: 
            fh.write("%s\t%s\n" % ( child, self.childToParent[child] ) )
        fh.close()




    def init_extractRandomlyFromTree(self, treeFile, maxLevels=7): 
        """Create a taxonomy that is consistent with a given tree. The
        maximum depth of a leave in this taxonomy is <maxLevel>. Also
        note, that the taxonomy will consider the given rooting for
        the tree. If the tree is unrooted the newick-string
        representation is used as a basis."""

        tree = Tree(treeFile)
        idfact = IdFactory( tree.get_leaf_names(), 0)
        self.root = idfact.produceId()
        self.__traversal(self.root, tree, maxLevels, 0, idfact)

        # populate the other dict 
        for child in  self.childToParent.keys():
            parent = self.childToParent[child]            
            if self.parentToChildren.has_key(parent): 
                self.parentToChildren[parent].append(child)
            else : 
                self.parentToChildren[parent] = [child]
        self.dirty = True
        self.cleanup()        


    def __traversal(self, currentParent, tree, maxLevels, currentDepth, idFact): 
        if(tree.is_leaf()): 
            pass 
        else : 
            if(random.random() < 0.5 and currentDepth < maxLevels): 
                newId = idFact.produceId()
                self.childToParent[newId] = currentParent
                currentParent = newId
                for leave in tree.get_leaf_names(): 
                    self.childToParent[leave] = newId
            for elem in tree.get_children(): 
                self.__traversal(currentParent, elem, maxLevels, currentDepth + 1, idFact)


    def cleanup(self): 
        " To be called when dirty. "
        if(self.dirty): 
            self.reduceTaxonomy(self.getLeaves())
        self.dirty = False


    def reduceTaxonomy(self, relevantLeaves): 
        " important: reduces the taxonomy (in place) to a set of leaves "
        self.__rebuildWithRelevant(set(relevantLeaves))
        self.__checkConsistency()


    def getNthBipartition(self, leave, n): 
        """
        Gets all leaves in the bipartiion that is obtained by going up
        <n>-levels from <leave>.
        """
        result = set()
        leaves = self.getLeaves()
        innerNode = leave

        assert(self.childToParent.has_key(innerNode))
        # go up 
        for i in xrange(0, n): 
            if innerNode == self.root: 
                break
            else : 
                innerNode = self.childToParent[innerNode]

        return self.__getBipHelper(innerNode, leaves, result)


    def mislabel(self, leave, numRankUp, numRankDown): 
        """mislabel taxon <leave>. Put this taxon <numRankUp> up and then <numRankDown> 

        There is a chance, that taxa will not get mislabelled that
        way. You should check on equality later with origTax ==
        mislabelTax, where origTax is the taxonomy before and
        mislabelTax is a deep copy on which the mislabel method has
        been called.
        """ 
        self.dirty = True
        
        leaves = self.getLeaves()            
        assert(leave in leaves)
        
        # go n steps up  
        previous = leave
        currentNode = leave 
        firstParent = self.childToParent[leave]
        for i in range(0,numRankUp): 
            if currentNode == self.root: 
                break
            previous = currentNode
            currentNode = self.childToParent[currentNode] 
        
        # remove old relationship 
        del self.childToParent[leave] 
        self.parentToChildren[firstParent].remove(leave)

        # go m steps down 
        for i in range(0, numRankDown): 
            if(currentNode in leaves ):
                currentNode = self.childToParent[currentNode]
                break
            else : 
                currentNode = random.choice(list(set(self.parentToChildren[currentNode]) - set(previous)))
        if(currentNode in leaves): 
            currentNode = self.childToParent[currentNode] 
        assert(currentNode != leave)
        assert(currentNode not in leaves)
        self.childToParent[leave] = currentNode
        self.parentToChildren[currentNode].append(leave)
        self.cleanup()


    def __eq__(self,other): 
        """
        Test equality with == 
        """

        assert(type(other) == Taxonomy)
        
        if self.getMaxLevel() != other.getMaxLevel(): 
            return  False

        assert(self.getMaxLevel() == other.getMaxLevel())
        leaves = self.getLeaves()
        md = self.getMaxLevel()

        isEqual = True
        for leave in leaves: 
            for level in range(1, md) : 
                isEqual = isEqual and ( set(self.getNthBipartition(leave, level)) == set(other.getNthBipartition(leave, level)) )         
        return isEqual


############
# GETTERS  #
############
    def getBipartitionComplement(self, bipartition): 
        " Gets the complement of the taxa in the set <bipartition> "
        return self.getLeaves() - set(bipartition)
    
    

    def getDepth(self, node): 
        " Get the depth of a node "
        if(node == self.root): 
            return 0
        else : 
            return self.getDepth(self.childToParent[node]) + 1 
        

    def getMaxLevel(self): 
        " Get the number of levels in the taxonomy "
        return max(map(lambda x : self.getDepth(x) , self.getLeaves()))


    def getRoot(self): 
        return self.root

    
    def getInnerNodes(self): 
        return set(self.parentToChildren.keys())

    def getLeaves(self): 
        parents = set(self.parentToChildren.keys())
        children = set(self.childToParent.keys())
        return children - parents



######################
## PRINTING METHODS  #
######################
    def getNewickString(self): 
        return self.__newickHelper(self.root) + ';'

    def getChildToParentString(self): 
        result = "ROOT:"  + self.root + "\n"
        isFirst = True
        for child in self.childToParent.keys(): 
            if isFirst: 
                isFirst = False
            else :
                result += "\n"
            result += child + "\t" + self.childToParent[child]
        return result  

    def getParentToChildrenString(self): 
        result = "ROOT:" + self.root + "\n"
        isFirst = True
        for  parent in self.parentToChildren.keys(): 
            if isFirst : 
                isFirst = False
            else: 
                result += "\n"
            result += parent + "\t" + ",".join(self.parentToChildren[parent]) 
        return result


############
# PRIVATE  #
############

    def __handleMonofurcations(self, node, leaves): 
        if node in leaves:  
            return node 
        elif len(self.parentToChildren[node] ) > 1 : 
            newChildren =  []
            for child in self.parentToChildren[node]: 
                newChild = self.__handleMonofurcations(child, leaves)
                newChildren.append(newChild)
                if(newChild != child): 
                    del self.childToParent[child]
                    self.childToParent[newChild] = node 
            self.parentToChildren[node] = newChildren
            return node 
        else:
            assert(len(self.parentToChildren[node] ) == 1 ) 
            currentChild = self.parentToChildren[node][0]

            # remove self 
            del self.parentToChildren[node]
            del self.childToParent[currentChild] 

            return self.__handleMonofurcations(currentChild, leaves)



    def __rebuildWithRelevant(self, relevantLeaves): 
        # update childToParent
        newChildToParent = {}
        for leave in relevantLeaves:
            currentNode = leave
            assert(currentNode != self.root)
            while(currentNode != self.root): 
                newChildToParent[currentNode] = self.childToParent[currentNode]
                currentNode = self.childToParent[currentNode]       
        self.childToParent = newChildToParent
        
        #  update parentToChildren
        for (parent, children) in self.parentToChildren.items():
            newChildren = []
            for child in children: 
                if self.childToParent.has_key(child): 
                    newChildren.append(child)
            if(len(newChildren) == 0 ): 
                del(self.parentToChildren[parent])
            else : 
                self.parentToChildren[parent] = newChildren

        # check from root on, if there are mono-furcations in the taxonomy 
        self.root = self.__handleMonofurcations(self.root, self.getLeaves())



    def __newickHelper(self, node): 
        if(self.parentToChildren.has_key(node)): 
            return "(" + ",".join(map(lambda x : self.__newickHelper( x), self.parentToChildren[node])) + ")"
        else : 
            return str(node)


        
    def __checkConsistency(self): 
        num = self.__getNumEdgesForCheck(self.root, self.getLeaves())
        assert(len(self.childToParent.items()) == num - 1 ) 


    def __getNumEdgesForCheck(self, node, leaves): 
        if(node in leaves): 
            return 1
        else : 
            result = 1 
        for child in self.parentToChildren[node]: 
            result +=  self.__getNumEdgesForCheck(child, leaves)
        return result 

            
    def __str__(self): 
        return self.getParentToChildrenString()         

    def __getBipHelper(self,innerNode, leaves, result): 
        if innerNode in leaves: 
            result.add(innerNode) 
        else : 
            for elem in self.parentToChildren[innerNode]: 
                self.__getBipHelper(elem, leaves, result) 
        return result


################
# END


if __name__ == "__main__":
    if len(sys.argv) != 2: 
        print "usage: ./script <file>"
        sys.exit()
        
    tax = Taxonomy()
    tax.init_extractRandomlyFromTree(sys.argv[1] )
    print tax.getNewickString()

    # tax = Taxonomy()
    # tax.init_parseTaxFile("test.tax")

    # tax2 = Taxonomy()
    # tax2.init_parseTaxFile("test2.tax")
    
    # if tax == tax2: 
    #     print "tax are same"
    # else : 
    #     print "tax are different"

