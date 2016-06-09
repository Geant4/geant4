# $Id: tiaraApplication.py,v 1.3 2003/06/16 17:06:44 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-05-02 $
# -------------------------------------------------------------------
#
import string

import CLHEP
import G4Kernel
import Tiara

import tiaraSpecifications
import tallyData



class TiaraApplet(object):
    tiaraSim = None
    def __init__(self, tiaraSpecs, tSim = None, useLizard = True):
        if useLizard:
            import myLiz
            self.tree = myLiz.tf.create()
            self.hf = myLiz.af.createHistogramFactory (self.tree)
        else:
            self.tree = None
            self.hf = None
            
        if (not TiaraApplet.tiaraSim):
            if not tSim:
                print "Error: TiaraApplet: argument is empty, and no tiaraSim exists!"
            else:
                TiaraApplet.tiaraSim = tSim
        self.tiaraSpecs = tiaraSpecs
        self.tiaraHall = Tiara.TiaraGeometry(self.tiaraSpecs.\
                                             materials)
        self.cellScorerStore = Tiara.TiaraCellScorerStore()
        self.sampler = G4Kernel.G4MassGeometrySampler("neutron")
        self.eventAction = None
        self.scorer = None
        self.primGen = None
        self.scSrc = None
        self.physicsList = None
        self.nameExt = ""
        self.scoreDets = []
        self.scoreDetectorCreator = None
        self.noComponents = 0
        self.createdDPS = 0
        self.tiaraSim.SetGeometry(self.tiaraHall)

    def setScoreDetectorCreator(self, creator):
        self.scoreDetectorCreator=creator
        
    def visMode(self):
        if self.eventAction:
            print "TiaraApplet.visMode(): event action exists already"
        else:
            self.eventAction = Tiara.TiaraVisEventAction()
            self.tiaraSim.AddVisRunAction()
            self.tiaraSim.AddTiaraEventAction(self.eventAction)
        return

    def timedMode(self, time, randomNumFileName = ""):
        if self.eventAction:
            print "TiaraApplet.timedMode(): event action exists already"
        else:
            self.eventAction = Tiara.TiaraTimedEventAction(time)
            if randomNumFileName:
                self.eventAction.SetRnadomNumFilename(randomNumFileName)
            self.tiaraSim.AddTiaraEventAction(self.eventAction)
        return
    

    def specifyPhysicsList(self, pl, particleCut):
        self.particleCut = particleCut
        self.physicsList = pl
    
    def buildColimator(self):
        posColi = \
                CLHEP.Hep3Vector(0,0, self.tiaraSpecs.
                                 dimensions.targetPosZ + 
                                 self.tiaraSpecs.
                                 dimensions.distTargetExperiment +
                                 self.tiaraSpecs.\
                                 experiment.colWidth/2)
        
        logCol = self.tiaraHall.BuildCollimator(self.tiaraSpecs.\
                                                experiment.colWidth,
                                                "iron",
                                                "air")
        self.tiaraHall.PlaceExpComponent(posColi,
                                         logCol,
                                         "colimator")
        return

    def buildShield(self):
        log = self.tiaraHall.BuildShield(self.tiaraSpecs.\
                                         experiment.shieldWidth,
                                         self.tiaraSpecs.\
                                         experiment.shieldMaterial)
        posShield = CLHEP.Hep3Vector(0,0, self.tiaraSpecs.
                                     dimensions.targetPosZ + 
                                     self.tiaraSpecs.
                                     dimensions.
                                     distTargetExperiment + 
                                     self.tiaraSpecs.experiment.
                                     colWidth +
                                     self.tiaraSpecs.experiment.\
                                     shieldWidth/2)
        self.tiaraHall.PlaceExpComponent(posShield, log, "shield")
        return

    def buildGeometry(self):
        self.tiaraHall.BuildGeometry(self.tiaraSpecs.dimensions)
        if self.noComponents != 1:
            self.tiaraHall.CreateComponents()
        if self.tiaraSpecs.experiment.colWidth > 0:
            self.buildColimator()
        self.buildShield()

    def createCellScorers(self):
        self.cellScorer = []
        tally = Tiara.TiaraTally()
        tally.setBinEdges(tiaraSpecifications.\
                          tallyBinEdges[self.tiaraSpecs.\
                                        experiment.energy])
        for det in self.scoreDets:
            det.scorer = None
            if self.hf:
                det.scorer = \
                           Tiara.\
                           TiaraCellScorer(self.hf,
                                           det.name,
                                           self.tiaraSpecs.
                                           experiment.binEdgesScinti,
                                           self.tiaraSpecs.
                                           experiment.binEdgesBonner,
                                           tally)
            else:
                det.scorer = \
                           Tiara.\
                           TiaraCellScorer(det.name,
                                           tally)
            self.cellScorer.append(det.scorer)
            
                                   
    def buildDetectors(self):
        tally = Tiara.TiaraTally()
        tally.setBinEdges(tiaraSpecifications.\
                          tallyBinEdges[self.tiaraSpecs.
                                        experiment.energy])

        physSrcDet = self.tiaraHall.AddSourceDetector()
        self.scSrc = None
        if self.hf:
            self.scSrc = Tiara.\
                         TiaraCellScorer(self.hf,
                                         "source_detector",
                                         self.tiaraSpecs.
                                         experiment.binEdgesScinti,
                                         self.tiaraSpecs.
                                         experiment.binEdgesBonner,
                                         tally)
        else:
            self.scSrc = Tiara.\
                         TiaraCellScorer("source_detector",
                                         tally)
            

        self.scoreDets = self.scoreDetectorCreator.\
                         createScoreDetectors(self.tiaraHall)

        self.createCellScorers()

        self.cellScorerStore.AddTiaraCellScorer(self.scSrc,
                                      G4Kernel.G4GeometryCell(physSrcDet,
                                                                  0))
        
        for det in self.scoreDets:
            self.cellScorerStore.AddTiaraCellScorer(det.scorer,
                                                    G4Kernel.\
                                                    G4GeometryCell(det.
                                                                   phys,
                                                                   0))

        self.scorer = G4Kernel.G4CellStoreScorer(self.cellScorerStore.
                                                 GetG4VCellScorerStore())

        self.eventAction.SetScorerStore(self.cellScorerStore)
            
        return

    def setPrimaryGenerator(self, primGen):
        self.primGen = primGen
        self.tiaraSim.SetPrimaryGenerator(self.primGen)
        return

    def fillShelve(self, shelveDB):
        t = self.primGen.GetTally()
        shelveDB["generatorTally"] = tallyData.createTallyDat(t)
        t = self.scSrc.GetTally()
        shelveDB["sourceDetectorTally"] = tallyData.createTallyDat(t)
        for det in self.scoreDets:
            t = det.scorer.GetTally()
            shelveDB[det.name + "Tally"] = tallyData.createTallyDat(t)
        shelveDB["runTime"] = self.eventAction.GetTotalProcessedTime()
        
    def config(self):

        self.tiaraSim.SetPhysicsList(self.physicsList)

        self.buildDetectors()
    

        self.tiaraSim.initialize()
        if len(self.particleCut) > 0:
            for p in self.particleCut:
                self.tiaraSim.AddParticleCut(p, self.particleCut[p]);
        
        self.sampler.PrepareScoring(self.scorer)
        self.sampler.Configure()



