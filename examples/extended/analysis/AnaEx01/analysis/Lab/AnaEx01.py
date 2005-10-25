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
# Example of a Python script using the AnaEx01.aida file.
#
#/////////////////////////////////////////////////////////////////////////////

#
#  To be executed from the Python shell :
#   OS> <setup Lab>
#   OS> python
#  (>>> import AIDA)
#   >>> import AnaEx01
#

from AIDA import AIDA_createAnalysisFactory

aida = AIDA_createAnalysisFactory()

# Get some AIDA factories :
plotterFactory = aida.createPlotterFactory()
treeFactory = aida.createTreeFactory()
memoryTree = treeFactory.createDefault()
histogramFactory = aida.createHistogramFactory(memoryTree)

# Attach the current plotter :
plotter = plotterFactory.create('')

plotter.createRegions(1,2,0)
plotter.setParameter('title','AnaEx01 analysis')

tree = treeFactory.create('AnaEx01.aida','xml',1,0)
tree.ls()

#/////////////////////////////////////////////////////////////////////////////
# In first region, get and plot the EAbs histo :
#/////////////////////////////////////////////////////////////////////////////

tree.cd('histograms')
tree.ls()

# From an AIDA tree we retreive IManagedObjects, we have to cast.
#FIXME : IManagedObject_to_IHistogram1D is not AIDA, it is OSC specific.
from AIDA import IManagedObject_to_IHistogram1D
EAbs = IManagedObject_to_IHistogram1D(tree.find('EAbs'))
if EAbs == None:
  print 'EAbs not found or is not an AIDA::IHistogram1D.'
else :
  # Plot the histo :
  region = plotter.currentRegion()
  region.plot(EAbs)

#/////////////////////////////////////////////////////////////////////////////
# In second region plot the EAbs histo built from the tuple :
#/////////////////////////////////////////////////////////////////////////////

tree.cd('..')
tree.cd('tuples')
tree.ls()

tupleFactory = aida.createTupleFactory(tree)
# Get a tuple from the tree :
#tuple = tupleFactory.create('AnaEx01','AnaEx01','')
from AIDA import IManagedObject_to_ITuple
tuple = IManagedObject_to_ITuple(tree.find('AnaEx01'))
if tuple == None:
  print 'Tuple AnaEx01 not found or is not an AIDA::ITuple.'

# Create an histo in the default memory tree :
tuple_EAbs = histogramFactory.createHistogram1D('tuple_EAbs','AnaEx01/EAbs',100,0,100)

# Create an Evaluator and a Filter object :
evaluator = tupleFactory.createEvaluator('EAbs')
evaluator.initialize(tuple)
filter = tupleFactory.createFilter('')
filter.initialize(tuple)

# Project tuple AnaEx01/EAbs column on the tuple_EAbs histo :
tuple.project(tuple_EAbs,evaluator,filter)

plotter.next()
region = plotter.currentRegion()
region.plot(tuple_EAbs)

plotter.show()
plotter.interact()

del evaluator
del filter

print 'The two histos must give the same information !'

del tupleFactory

del tree
del plotter
del region
del plotterFactory
del treeFactory
del histogramFactory
