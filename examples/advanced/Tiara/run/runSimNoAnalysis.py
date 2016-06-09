#!/usr/bin/env python
#
# $Id: runSimNoAnalysis.py,v 1.6 2005/03/17 19:48:27 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-07-01 $
# -------------------------------------------------------------------


import CLHEP
import G4Kernel
import Tiara
import tiaraApplication
import tiaraGenerators
import tiaraDetectors
import tiaraSpecifications
import myUtils
import variableGeometry
import slabGeometry
import runSequence
import os

##########################################################################
# random number initialization
##########################################################################
Tiara.setRandomSeed(891011);
#Tiara.setRandomStatus("dTest/tiara-2003_5_27_20_5_39_pcgeant2/randomNumberFile_run00003");    










##########################################################################
# experiment and simulation specific data 
##########################################################################
particleCut = {"neutron" : 3 * CLHEP.MeV,
               "gamma"   : 1 * CLHEP.MeV,
               "proton"  : 1 * CLHEP.MeV,
               "deuteron": 1 * CLHEP.MeV,
               "triton"  : 1 * CLHEP.MeV,
               "alpha"   : 1 * CLHEP.MeV}

beamEnergy = 43
shieldWidth = 150 * CLHEP.cm

totalTime = 3 * myUtils.min
timeForOneRun = 1 * myUtils.min

# available physics lists: LHEP_LEAD, LHEP_PRECO_HP
# LHEP_PRECO, LHEP
##physList = Tiara.LHEP()
##physList = Tiara.LHEP_PRECO()
physList = Tiara.LHEP_LEAD()
##physList = Tiara.LHEP_PRECO_HP()

# specify the detectors
scoreDetectorCreator = tiaraDetectors.ThreeZylindricDetectors()
#scoreDetectorCreator = tiaraDetectors.DetectorSlab()


comment = ""



##########################################################################
# Create a Specification object of the configuaration data
##########################################################################
experiment = tiaraSpecifications.Experiment(beamEnergy,      
                                           particleCut["neutron"],
                                           particleCut,
                                           shieldWidth,
                                           "concrete")

tiaraSpecs = tiaraSpecifications.Specifications(Tiara.TiaraDimensions(),
                                                experiment,
                                                Tiara.TiaraMaterials())



##########################################################################
# definition of the importance geometry and a scorer 
##########################################################################
impGeo = variableGeometry.VariableImpSlabGeometry(tiaraSpecs)

impGeo.addCellImportance(width=15.0 * CLHEP.cm, faktor=1)
for i in range(9):
    impGeo.addCellImportance(width=15.0 * CLHEP.cm, faktor=2)

impGeo.construct()

# an alternative
#impGeo = slabGeometry.SlabedImportanceGeometry(tiaraSpecs,
#                                               10.0 * CLHEP.cm,
#                                               1.5)


impScorer = G4Kernel.G4Scorer()



##########################################################################
# Creation of a TiaraApplet to define the run mode, physics list,
# detector type and the primary generator
##########################################################################
tApp = tiaraApplication.TiaraApplet(tiaraSpecs = tiaraSpecs,
                                    tSim = Tiara.TiaraSim_GetTiaraSim(),
                                    usePI = False)


#tApp.visMode()
tApp.timedMode(timeForOneRun)


tApp.specifyPhysicsList(physList, particleCut)

# detectors
tApp.setScoreDetectorCreator(scoreDetectorCreator)

tApp.buildGeometry()

tiara_dir = os.environ["TIARA_BASE"]

#primGenBuilder = tiaraGenerators.\
#                 TiaraDPSEnergyGenerator(tiaraSpecs,
#                                         tiara_dir +
#                                         "/data/expDataConverted/dpsSource.xml")
#primGenBuilder = tiaraGenerators.TiaraPrimaryGenerator(tiaraSpecs)
primGenBuilder = tiaraGenerators.FixedEnergyPrimaryGenerator(tiaraSpecs)

tApp.setPrimaryGenerator(primGenBuilder.primGen)


tApp.noComponents = 0

## The following lines should be probably moved here after release Geant4-7.0
##physList = Tiara.LHEP()
##physList = Tiara.LHEP_PRECO()
##physList = Tiara.LHEP_PRECO_HP()
##physList = Tiara.LHEP_LEAD()
##tApp.specifyPhysicsList(physList, particleCut)
##tApp.setPrimaryGenerator(primGenBuilder.primGen)

tApp.config()




##########################################################################
# creating of the sampler
##########################################################################
parallelSampler = myUtils.createParallelSampler(impGeo,
                                                impScorer)




##########################################################################
# create a run config object, a run sequence and run the simulation
##########################################################################
rc = runSequence.RunConfig()

rc.basePath = "simData"
rc.tApp = tApp
rc.tiaraSpecs = tiaraSpecs
rc.impGeo = impGeo
rc.impScorer = impScorer
rc.totalTime = totalTime
rc.comment = comment

rs = runSequence.RunSequence(runConfig=rc, usePI = False)
rs.runNevents(100)
rs.runLoop()


##########################################################################
##########################################################################
