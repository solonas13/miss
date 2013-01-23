#! /usr/bin/python

import random
import sys 
from ete2 import *

if len(sys.argv) != 3 : 
    print "./script <tree> <numTax>"
    sys.exit()

tr = Tree(sys.argv[1])
numTax = int(sys.argv[2])
leaves = tr.get_leaf_names()
random.shuffle(leaves)
toPrune = set(leaves) - set(leaves[0:numTax])
tr.prune(toPrune)
print tr.write(format=2)
