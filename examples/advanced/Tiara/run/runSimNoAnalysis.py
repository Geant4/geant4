#!/usr/bin/env python2.2
#
# $Id: runSimNoAnalysis.py,v 1.3 2003/06/20 12:41:06 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-06-00 $
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

# available physics lists: TiaraPhysicsList, LHEP_LEAD_HP, LHEP_PRECO_HP
# CASCADE_HP LHEP_BIC  LHEP_BIC_BIC
physList = Tiara.LHEP_BIC_HP()

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
                                    useLizard = False)


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

rs = runSequence.RunSequence(runConfig=rc, useLizard = False)
rs.runNevents(100)
rs.runLoop()


##########################################################################
##########################################################################
