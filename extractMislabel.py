#!/usr/bin/python

import os
import random
import re
import sys

def label_exists ( tree, label ):
  x = re . search ( label + "[\,\)]", tree )
  if x == None:
    return False
  else:
    return True

def read_raxml_output ( filename ):
  s = open ( filename ) . read ()
  s = s . split ( '\n' )

  return ( s[1], s[3:] )    # return reference tree and placement strings

def normalise_reference_tree ( rt ):
  rt = re . sub ( r"^[^(]*","", rt )
  rt = re . sub ( r":[^,()]*", "", rt )
  rt = rt[:-3]

  return ( rt )

def get_node_of_label ( tree, label ):
  node = re . search ( label + ":" + "\d*\.\d*\{(\d+)\}", tree )
  return ( node . group ( 1 ) )

def get_label_of_node ( tree, node ):
  label = re . search ( "([\w]+):\d*\.\d*\{" + node + "\}", tree )
  if ( label != None ):
    return ( label . group ( 1 ) )
  print "seq Label not found for node " + node + " , internal branch"
  return False 

def get_max_lh_node ( placement, node ):
  #x = re . search ( "\[" + node + "(\,\s*\-?\d*\.?\d*){3}\,\s*(\-?\d*\.?\d*)\]", placement )
  #x = re . search ( "\[" + node + "\,\s*(\-?\d*\.?\d*)", placement )
  #return ( x . group ( 1 ) )
  x = re . search ( "\[" + node + "\,\s*(\-?\d*\.?\d*),\s*(\d+\.\d+)", placement )
  return  " LH " + x . group ( 1 ) + " Weight: " +  x . group ( 2 ) 

# Get the node with the maximum likelihood
def get_max_lh ( placement ):
  #x = re . search ( ":\[\[(\d+),\s*([^,]+)", placement )
  x = re . search ( ":\[\[(\d+),\s*([^,]+),\s*(\d+\.\d+)", placement )
  if ( x != None ):
    return ( x . group ( 1 ), x . group ( 2 ), x . group ( 3 ) )
  return False


def add_sibling ( tree, where, node ):
  ltree = re . findall ( r"[\w]+|[\(\)\,\;]", tree )
  i = ltree . index ( where )
  ltree . insert ( i + 1, ')' )
  ltree . insert ( i, ',' )
  ltree . insert ( i, node )
  ltree . insert ( i, '(' )
  return "" . join ( ltree )

def get_siblings ( tree, x ):
  p = re . search ( x + "[\,\)]", tree )
  if p == None:
    print "No such node found"
    return

  y = tree . partition ( p . group ( 0 ) )
  siblings = []

  if y[1] == y[2] == "":
    print "Separator not found"
    return
  
  pre  = re . findall ( r"[\w]+|[\(\)\,\;]", y[0] )
  post = re . findall ( r"[\w]+|[\(\)\,\;]", y[2] )

  level = 0
  for c in reversed ( pre ):
    if c == ')':
      level += 1
    elif c == '(':
      if level == 0:
        break;
      level -= 1
    elif c != ',':
      if level == 0:
        siblings . append ( c )

  if y[1][-1:] != ')':
    level = 0
    for c in post:
      if c == '(':
        level += 1
      elif c == ')':
        if level == 0:
          break;
        level -= 1
      elif c != ',':
        if level == 0:
          siblings . append ( c )
    
  return siblings

def get_inserted_label ( p ):
  label = re . search ( "\"n\":\[\"(.*)\"\]", p );
  return label . group ( 1 );
 

def main ( argv = None ):
  if argv is None: argv = sys . argv

  if len ( argv ) != 3:
    print "usage: ./tree.py <RAxML.out> <taxonomy.out>";
    sys . exit ( 0 )

  # get the reference tree and placements
  ( rt, placements ) = read_raxml_output ( argv[1] )
  ref_tree = normalise_reference_tree ( rt )
  print ref_tree

  # read and print the taxonomy tree
  tax_tree = open ( argv[2] ) . read ()
  print tax_tree . rstrip ( '\n' )

  for p in placements:
    if get_max_lh ( p ) == False:
      return
    (node, L1, W1 ) = get_max_lh ( p )

    print
    new_label = get_inserted_label ( p )
    print "\nnew label (taxon to be added) " + new_label
    label = get_label_of_node ( rt, node )
    if label != False:
      print "ML EPA-insertion node " + node + " has label " + label
      tree = add_sibling ( ref_tree, label, new_label )
      #print tree
    print "Max Likelihood placement at node " + node +  " LH " + L1 + ", weight " + W1

    siblings = get_siblings( tax_tree, new_label )
    print siblings
    for label in siblings:
      if label_exists ( ref_tree, label ) == False:
        continue
      tree = add_sibling ( ref_tree, label, new_label )
      #print tree
      node = get_node_of_label ( rt, label )
      print "Taxonomy aware placement at " + label + get_max_lh_node ( p, node )

if __name__ == "__main__":
  sys . exit ( main ( ) )


