#!/usr/bin/env python2.2
import Tiara
import G4Kernel
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
particleCut = {"neutron" : 3 * G4Kernel.MeV,
               "gamma"   : 1 * G4Kernel.MeV,
               "proton"  : 1 * G4Kernel.MeV,
               "deuteron": 1 * G4Kernel.MeV,
               "triton"  : 1 * G4Kernel.MeV,
               "alpha"   : 1 * G4Kernel.MeV}

beamEnergy = 43
shieldWidth = 150 * G4Kernel.cm

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
# definition of the importance geometry importance values and a scorer 
##########################################################################
impGeo = variableGeometry.VariableImpSlabGeometry(tiaraSpecs)

impGeo.addCellImportance(width=15.0 * G4Kernel.cm, faktor=1)
for i in range(9):
    # to run unbiased set: faktor=1     in the next line
    impGeo.addCellImportance(width=15.0 * G4Kernel.cm, faktor=2)

impGeo.construct()

# an alternative
#impGeo = slabGeometry.SlabedImportanceGeometry(tiaraSpecs,
#                                               10.0 * G4Kernel.cm,
#                                               1.5)


impScorer = G4Kernel.G4Scorer()







##########################################################################
# Creation of a TiaraApplet to define the run mode, physics list,
# detector type and the primary generator
##########################################################################
tApp = tiaraApplication.TiaraApplet(tiaraSpecs,
                                    Tiara.TiaraSim_GetTiaraSim())



#tApp.visMode()
tApp.timedMode(timeForOneRun)


tApp.specifyPhysicsList(physList, particleCut)

# detectors
tApp.setScoreDetectorCreator(scoreDetectorCreator)

tApp.buildGeometry()

tiara_dir = os.environ["TIARA_BASE"]

primGenBuilder = tiaraGenerators.\
                 TiaraDPSEnergyGenerator(tiaraSpecs,
                                         tiara_dir +
                                         "/data/expDataConverted/dpsSource.xml")
#primGenBuilder = tiaraGenerators.TiaraPrimaryGenerator(tiaraSpecs)
#primGenBuilder = tiaraGenerators.FixedEnergyPrimaryGenerator(tiaraSpecs)

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

rs = runSequence.RunSequence(rc)
rs.runNevents(100)
rs.runLoop()


##########################################################################
##########################################################################
