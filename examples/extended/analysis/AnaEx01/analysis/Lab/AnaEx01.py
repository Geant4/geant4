#
# ********************************************************************
# * DISCLAIMER                                                       *
# *                                                                  *
# * The following disclaimer summarizes all the specific disclaimers *
# * of contributors to this software. The specific disclaimers,which *
# * govern, are listed with their locations in:                      *
# *   http://cern.ch/geant4/license                                  *
# *                                                                  *
# * Neither the authors of this software system, nor their employing *
# * institutes,nor the agencies providing financial support for this *
# * work  make  any representation or  warranty, express or implied, *
# * regarding  this  software system or assume any liability for its *
# * use.                                                             *
# *                                                                  *
# * This  code  implementation is the  intellectual property  of the *
# * GEANT4 collaboration.                                            *
# * By copying,  distributing  or modifying the Program (or any work *
# * based  on  the Program)  you indicate  your  acceptance of  this *
# * statement, and all its terms.                                    *
# ********************************************************************
#
#/////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////
#
# Example of a Python script using the AnaEx01.root file.
#
#/////////////////////////////////////////////////////////////////////////////

# Hook the current plotter :
plotter = Plotter('')

plotter.createRegions(1,2)
plotter.setParameter('pageTitle','AnaEx01 analysis')

rootTree = RootTree('AnaEx01.root','READ')
rootTree.ls()

#/////////////////////////////////////////////////////////////////////////////
# In first region, get and plot the EAbs histo :
#/////////////////////////////////////////////////////////////////////////////

rootTree.cd('histograms')
rootTree.ls()

#  Get some histograms with the H1Get 
# builtin Python procedure (defined in Lab/scripts/Python/Lab_init.py).

EAbs = H1D_get(rootTree,'EAbs')

# Plot the histo :
plotter.plot(EAbs)

#/////////////////////////////////////////////////////////////////////////////
# In second region plot the EAbs histo built from the tuple :
#/////////////////////////////////////////////////////////////////////////////

rootTree.cd('..')
rootTree.cd('tuples')
rootTree.ls()

tuple = Tuple("AnaEx01","AnaEx01","",rootTree)

# Create an histo in the default memory tree :
tuple_EAbs = H1D('tuple_EAbs','AnaEx01/EAbs',100,0,100)

# Create an Evaluator and a Filter object :
evaluator = Evaluator(tuple,'EAbs')
filter = Filter(tuple,'')

# Project tuple AnaEx01/EAbs column on the tuple_EAbs histo :
tuple.project(tuple_EAbs,evaluator,filter)

plotter.next()
plotter.plot(tuple_EAbs)
plotter.setParameter('histogramContext','color red modeling solid')

# Evaluators and Filters are not managed, delete them explicitly :
delete_IEvaluator(evaluator)
delete_IFilter(filter)

echo('The two histos must give the same information !')
