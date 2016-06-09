#!/usr/bin/env python
#
# $Id: runSim.py,v 1.12 2005/12/15 16:23:11 ahoward Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-08-01 $
# -------------------------------------------------------------------

# importing python libraries
import os

# importing wrapper modules (see the source/*Wrapper directories)
# created by swig (from the source/*Wrapper/*.i files).
import CLHEP
import G4Kernel
import Tiara

# importing python modules specific to this example 
# (see source/py_modules)
import tiaraApplication
import tiaraGenerators
import tiaraDetectors
import tiaraSpecifications
import myUtils
import variableGeometry
import slabGeometry
import runSequence

##########################################################################
# random number initialization
##########################################################################
Tiara.setRandomSeed(891011);
#Tiara.setRandomStatus("dTest/tiara-2003_5_27_20_5_39_pcgeant2/randomNumberFile_run00003");    





##########################################################################
# experiment and simulation specific data 
##########################################################################
# create a list of particles and give a lower energy cut
particleCut = {"neutron" : 3 * CLHEP.MeV,
               "gamma"   : 1 * CLHEP.MeV,
               "proton"  : 1 * CLHEP.MeV,
               "deuteron": 1 * CLHEP.MeV,
               "triton"  : 1 * CLHEP.MeV,
               "alpha"   : 1 * CLHEP.MeV}

# specify if the source neutrons shoul be form the 43 or 68 MeV protons
beamEnergy = 68 * CLHEP.MeV
# specify the shieldwidth [25, 50, 100, 150, (200 for 68 MeV case only)]
shieldWidth = 100 * CLHEP.cm

# set the total running time and the time for one run
# followed by a printout
# if the total running time, giving 0 will start a visualization
# and timeForOneRun has no meaning
totalTime = 5 * myUtils.min
timeForOneRun = 2 * myUtils.min
# if running the visualization mode by setting totalTime = 0,
# a Geant4 session will be started. Start the visualization
# by e.g.:
# Idle> control/execute vis.mac
# Idle> /run/beamOn 10


# available physics lists: LHEP_LEAD, LHEP_PRECO_HP
# LHEP_BIC  LHEP_BIC_HP  LHEP_PRECO  LHEP
#3physList = Tiara.LHEP()
##physList = Tiara.LHEP_PRECO()
##physList = Tiara.LHEP_PRECO_HP()
physList = Tiara.LHEP_LEAD()


# specify the detectors
scoreDetectorCreator = tiaraDetectors.ThreeZylindricDetectors()
#scoreDetectorCreator = tiaraDetectors.DetectorSlab()


comment = ""




##########################################################################
# Create a Specification object of the configuaration data
##########################################################################

# this should be fine for most settings
experiment = tiaraSpecifications.Experiment(beamEnergy,      
                                           particleCut["neutron"],
                                           particleCut,
                                           shieldWidth,
                                           "concrete")

tiaraSpecs = tiaraSpecifications.Specifications(Tiara.TiaraDimensions(),
                                                experiment,
                                                Tiara.TiaraMaterials())




##########################################################################
# definition of the importance geometry importance values and a scorer 
##########################################################################
# create a parallel geometry.
impGeo = variableGeometry.VariableImpSlabGeometry(tiaraSpecs)


# introduce "importance cells" into the shielding region
# In this case the width of all the shields in the shielding region
# are equal and the importance staring from one doubles
# from cell to cell in the beam direction.
impGeo.addCellImportance(width=20.0 * CLHEP.cm, faktor=1)
for i in range(4):
    # to run unbiased set: faktor=1     in the next line
    impGeo.addCellImportance(width=20.0 * CLHEP.cm, faktor=2)

impGeo.construct()

# Make sure the widths given in the above way add up to the
# shield widths e.g. 25, 50, ... cm.
# The importance geometries take care to give the importance value 1
# to the volume before the shield and to assign the same importance
# as the last cell in the shield region to the volume behind the shield.

# an alternative
#impGeo = slabGeometry.SlabedImportanceGeometry(tiaraSpecs,
#                                               10.0 * CLHEP.cm,
#                                               1.5)


impScorer = G4Kernel.G4Scorer()



##########################################################################
# Creation of a TiaraApplet to define the run mode, physics list,
# detector type and the primary generator
##########################################################################

# this and the remaining part should be fine for most settings
tApp = tiaraApplication.TiaraApplet(tiaraSpecs,
                                    Tiara.TiaraSim_GetTiaraSim())

if totalTime == 0:
    tApp.visMode()
else:
    tApp.timedMode(timeForOneRun)


tApp.specifyPhysicsList(physList, particleCut)

# detectors
tApp.setScoreDetectorCreator(scoreDetectorCreator)
tApp.buildGeometry()
tiara_dir = os.environ["TIARA_BASE"]

tApp.setPhysics()

primGenBuilder = tiaraGenerators.\
                 TiaraDPSEnergyGenerator(tiaraSpecs,
                                         tiara_dir +
                                         "/data/expDataConverted/dpsSource.xml")
#primGenBuilder = tiaraGenerators.TiaraPrimaryGenerator(tiaraSpecs)
#primGenBuilder = tiaraGenerators.FixedEnergyPrimaryGenerator(tiaraSpecs)

tApp.setPrimaryGenerator(primGenBuilder.primGen)

tApp.noComponents = 0
tApp.config()

## The following lines should be probably moved here after release Geant4-7.0
##physList = Tiara.LHEP()
##physList = Tiara.LHEP_PRECO()
##physList = Tiara.LHEP_PRECO_HP()
##physList = Tiara.LHEP_LEAD()
##tApp.specifyPhysicsList(physList, particleCut)
##tApp.setPrimaryGenerator(primGenBuilder.primGen)



##########################################################################
# creating of the sampler
##########################################################################
parallelSampler = myUtils.createParallelSampler(impGeo,
                                                impScorer)



##########################################################################
# create a run config object, a run sequence and run the simulation
##########################################################################
if totalTime > 0:
    rc = runSequence.RunConfig()

    rc.basePath = "simData"
    rc.tApp = tApp
    rc.tiaraSpecs = tiaraSpecs
    rc.impGeo = impGeo
    rc.impScorer = impScorer
    rc.totalTime = totalTime
    rc.comment = comment

    rs = runSequence.RunSequence(rc)
    rs.runNevents(100)
    rs.runLoop()
else:
    tApp.tiaraSim.startSession()


##########################################################################
##########################################################################
