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

from Lab import *

# Get some AIDA factories :
plotterFactory = aida.createPlotterFactory()
treeFactory = aida.createTreeFactory()
histogramFactory = aida.createHistogramFactory(memoryTree)
functionFactory = aida.createFunctionFactory(memoryTree)

# Attach the current plotter :
plotter = plotterFactory.create('')

plotter.createRegions(1,2,0)
plotter.setParameter('title','AnaEx01 analysis')

# If already loaded close (not delete) the file :
fileName='AnaEx01.root'
session.destroyManager(fileName)

rioTree = treeFactory.create(fileName,'ROOT',1,0)
rioTree.ls()

#/////////////////////////////////////////////////////////////////////////////
# In first region, get and plot the EAbs histo :
#/////////////////////////////////////////////////////////////////////////////

rioTree.cd('histograms')
rioTree.ls()

#  Get some histograms with the H1Get 
# builtin Python procedure (defined in Lab/scripts/Python/Lab_init.py).

EAbs = IManagedObject_to_IHistogram1D(rioTree.find('EAbs'))
if EAbs == None:
 print 'EAbs not found or is not an AIDA::IHistogram1D.'
else :
 # Plot the histo :
 region = plotter.currentRegion()
 region.plot(EAbs)

#/////////////////////////////////////////////////////////////////////////////
# In second region plot the EAbs histo built from the tuple :
#/////////////////////////////////////////////////////////////////////////////

rioTree.cd('..')
rioTree.cd('tuples')
rioTree.ls()

rioTupleFactory = aida.createTupleFactory(rioTree)
# Get a tuple from the rioTree :
tuple = rioTupleFactory.create('AnaEx01','AnaEx01','')

# Create an histo in the default memory tree :
tuple_EAbs = histogramFactory.createHistogram1D('tuple_EAbs','AnaEx01/EAbs',100,0,100)

# Create an Evaluator and a Filter object :
evaluator = rioTupleFactory.createEvaluator('EAbs')
evaluator.initialize(tuple)
filter = rioTupleFactory.createFilter('')
filter.initialize(tuple)

# Project tuple AnaEx01/EAbs column on the tuple_EAbs histo :
tuple.project(tuple_EAbs,evaluator,filter)

plotter.next()
region = plotter.currentRegion()
region.plot(tuple_EAbs)
region.setParameter('histogramContext','color red modeling solid')

del evaluator
del filter

print 'The two histos must give the same information !'

del rioTupleFactory

del rioTree
del plotter
del region
del plotterFactory
del treeFactory
del histogramFactory
del functionFactory
