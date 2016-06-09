#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   gtest04
#   - test for command tree and command information
# ==================================================================
from Geant4 import *

def DumpTree(atree):
  print "@@", atree.GetPathName(), "::", atree.GetTitle()
  ntree= atree.GetTreeEntry()
  ncommand= atree.GetCommandEntry()

  for i in range(1, ncommand+1):
    icommand= atree.GetCommand(i)
    print "  **", icommand.GetCommandPath()
    print "     ", icommand.GetTitle()
    x= icommand.GetStateList()

    nparameter= icommand.GetParameterEntries()
    for j in range(0, nparameter):
      iparam= icommand.GetParameter(j)
      print "    +", iparam.GetParameterName(), iparam.GetParameterType()
      
  for i in range(1, ntree+1):
    itree= atree.GetTree(i)
    DumpTree(itree)
    
# ==================================================================
# main
# ==================================================================
root_tree= gUImanager.GetTree()

DumpTree(root_tree)

