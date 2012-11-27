#!/usr/bin/python

import os
import random
import re
import sys

def parse_raxml_output ( filename ):
  f = open ( filename, 'r' )
  s = f . read ()

  s = s . split ( '\n' );
  rt = s[1]
  pl = s[3]

  return [rt, pl]

def construct_reference_tree ( rt ):
  # Print the reference tree RT
  rt = re . sub ( r"^[^(]*","", rt )
  rt = re . sub ( r":[^,()]*", "", rt )
  rt = rt[:-3]

  return ( rt )

def maximum_likelihood ( p ):
  x = re . search ( ":\[\[(\d+),\s*([^,]+)", p )

  return ( x . group ( 1 ), x . group ( 2 ) )

def find_inserted_label ( p ):
  label = re . search ( "\"n\":\[\"(.*)\"\]", p );
  return label . group ( 1 );

def compute_brackets_to_label ( rt, label ):
  tmp = re . search ( "(.*)" + label, rt );
  tmp = tmp . group ( 1 );
  print tmp;

  opening = len ( re . findall ( r"\(", tmp ) )
  closing = len ( re . findall ( r"\)", tmp ) )
  commas  = len ( re . findall ( r"\,", tmp ) )

  return ( opening, closing, commas )

def compute_brackets_to_node ( x, node ):
  tmp = re . search ( "(.*)\{" + str ( node ) + "\}", x );
  tmp = tmp . group ( 1 );
  opening = len ( re . findall ( r"\(", tmp ) )
  closing = len ( re . findall ( r"\)", tmp ) )
  commas  = len ( re . findall ( r"\,", tmp ) )
  return ( opening, closing, commas )

def insert_node ( label, ins_open, ins_close, commas, rt ):
 for i in range ( len ( rt ) ):
   if rt[i] == '(':
     ins_open -= 1
   elif rt[i] == ')':
     ins_close -= 1
   elif rt[i] == ',':
     commas -= 1
   if ins_open == 0 and ins_close == 0 and commas == 0:
     break
 

 i += 1;
 leaf = 0
 if rt[i] != '(' and rt[i] != ')' and rt[i] != ',':
   leaf = 1;
   open_bracket_pos = i;
 
 if leaf == 1:
   while rt[i] != ',' and rt[i] != ')':
     i += 1;
   rt = rt[:i] + ',' + label + ')' + rt[i:]
   rt = rt[:open_bracket_pos] + '(' + rt[open_bracket_pos:]
   return rt;
   #print rt
 else:
   rt = rt[:i] + ',' + label + ')' + rt[i:]

   stack = 0;
   i -= 1;

   stack = 1
   while stack != 0:
    i -= 1
    if rt[i] == '(':
      stack -= 1;
    elif rt[i] == ')':
      stack += 1
   rt = rt[:i] + '(' + rt[i:]
   return rt;
   #print rt
 
def main ( argv = None ):
  if argv is None: argv = sys . argv

  if len ( argv ) != 3:
    print " usage: ./tree.py <RAxML.out> <taxonomy.out>";
    sys . exit ( 0 )

  x  = parse_raxml_output ( argv[1] )
  rt = construct_reference_tree ( x[0] )
  (node,L1) = maximum_likelihood ( x[1] )
  label = find_inserted_label ( x[1] )

  print "RT = " + rt
  print "Label = " + label
  print "Node of max likelihood = " + str ( node )

  ( ins_open, ins_close, ins_comma ) = compute_brackets_to_node ( x[0], node )
  
  #print "Ins_open:  " + str ( ins_open )
  #print "Ins_close: " + str ( ins_close )
  #print "Comma: " + str ( ins_comma )

  t1 = insert_node ( label, ins_open, ins_close, ins_comma, rt )
  
  print "T1 = " + t1
  print "L1 = " + L1
  
  #( opening, closing, comma ) = compute_brackets_to_label ( rt, "Seq21" )
  #print "Opening brackets: " + str ( opening )
  #print "Closing brackets: " + str ( closing )



if __name__ == "__main__":
  sys . exit ( main ( ) )
