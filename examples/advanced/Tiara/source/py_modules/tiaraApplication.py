# $Id: tiaraApplication.py,v 1.7 2006/06/26 10:13:14 ahoward Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-08-02 $
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
    def __init__(self, tiaraSpecs, tSim = None, usePI = True):
        if usePI:
            import myPI
            self.tree = myPI.tf.create()
            self.hf = myPI.af.createHistogramFactory (self.tree)
            print self.hf
            print "hf"
            hf_ptr = str(self.hf._theObject)
            print hf_ptr
            print "hf_ptr before split"
            hf_ptr = string.split(str(self.hf._theObject))[-1]
            print hf_ptr
            print "hf_ptr after split" 
	    hf_ptr2 = hf_ptr.split('x')[1]
            for i in range(len(hf_ptr2),8) :
                hf_ptr2 = '0'+hf_ptr2
                
       	    self.hf.this='_' + hf_ptr2 + '_p_AIDA__IHistogramFactory'
            self.hf.thisown = 1
            print self.hf.this
            print "hf.this after split"
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

    def specifyPhysicsList(self, pl, particleCut):
        self.particleCut = particleCut
        self.physicsList = pl
        
    def setScoreDetectorCreator(self, creator):
        self.scoreDetectorCreator=creator
        
    def visMode(self):
        if self.eventAction:
            print "TiaraApplet.visMode(): event action exists already"
        else:
            self.eventAction = Tiara.TiaraVisEventAction()
# now later            self.tiaraSim.AddVisRunAction()
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

    def setPhysics(self):
        self.tiaraSim.SetPhysicsList(self.physicsList)
        return

    def visAdd(self):
        self.tiaraSim.AddVisRunAction()
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

        self.buildDetectors()
        self.tiaraSim.initialize()

        if len(self.particleCut) > 0:
            for p in self.particleCut:
                self.tiaraSim.AddParticleCut(p, self.particleCut[p]);
        
        self.sampler.PrepareScoring(self.scorer)
        self.sampler.Configure()

 
