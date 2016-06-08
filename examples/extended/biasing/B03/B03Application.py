# This file is used by G03RunApplication.py
#
# A very prototype like class which allows to do a given number
# of runs runs (nruns) of a given number of events (nevents)
# by calling nruns times BeamOn(nevents).
# The class fills a list after every call to BeamOn(nevents)
# with the number of tracks entered in a certain cell
# in the last run.



import math
import string

# import the G4 shadow classes. B03App in tern loads B03Appc.so
# containing the "hooks" to the c++ code 
from  B03App import *



class B03Application:
    def __init__(self, nevents):
        # this function prepares Geant4 for running by calling,
        # creates a parallel geometry, sets importance values
        # creates an parallel importance score sampler,
        # creates scores to access data,
        # and creates an empty histogram
        #(so by far to much for one function, but this is a prototype)
        

        # nevents is the number of events per run
        self.nevents = nevents
        # an array to store the messurments
        self.messurements = []

    def initializeApplication(self):
        if "base" in dir(self):
            print "already initialized"
            return
        # get B03AppBase singleton, in tern creating the mass detector and
        # the physics list 
        self.base = B03AppBase_GetB03AppBase ()


    def setupParallelGeometry(self):
        if "pgeom" in dir(self):
            print "geometry already created"
            return
        if not "base" in dir(self):
            print "call \"initializeApplication\" first"
            return
        # create the parallel geometry and get it's world volume
        self.pgeom = B03ImportanceDetectorConstruction()
        self.p_world = self.pgeom.GetWorldVolume()
        
        # create a list of G4GeometryCells related to the
        # physical volumes in the parallel detector
        # needed to assign importance values to the cells
        # and to create scorers for the cells
        self.gCell_list =[]
        for i in range(1,19):
            self.gCell_list.append(self.pgeom.GetGeometryCell(i))

        
        self.CellRest = self.pgeom.GetGeometryCell(19)

    def createIstoreAndScorer(self):
        if "istore" in dir(self):
            print "Istore and Scorer already created"
            return
        if not "pgeom" in dir(self):
            print "call \"setupParallelGeometry\" first"
            return
        
        # create an importance store  
        self.istore = G4IStore(self.p_world)

        # create a store for G4CellScorers
        self.cs_store = G4CellScorerStore()
        
        # fill the istore with importance values for the
        # G4GeometryCells
        # and create G4CellScorers for the G4GeometryCells 
        self.CellScorerList = []
        exp = 0
        for gCell in self.gCell_list:
            exp += 1
            importance = math.pow(2,exp)
            if exp < 18:
                self.istore.AddImportanceGeometryCell(importance,gCell)
            self.CellScorerList.append(self.cs_store.AddCellScorer(gCell))

        # set importance of world volume
        self.istore. \
        AddImportanceGeometryCell(1, G4GeometryCell(self.p_world,-1))
        imp_cell_17 = self.istore.GetImportance(self.pgeom.GetGeometryCell(17))
        
        self.istore. \
        AddImportanceGeometryCell(imp_cell_17,
                                  self.pgeom.GetGeometryCell(18))
                                  
        self.istore. \
        AddImportanceGeometryCell(imp_cell_17, self.CellRest)
            
        # create G4CellStoreScorer derived from G4VPScorer 
        self.scorer = G4CellStoreScorer(self.cs_store)

    
    def setupSampler(self):
        if "sampler" in dir(self):
            print "Sampler already created"
            return
        if not "istore" in dir(self):
            print "call \"createIstoreAndScorer\" first"
            return
        # create a sampler for neutrons
        self.sampler = G4ParallelGeometrySampler(self.p_world,
                                                 "neutron")

        self.sampler.PrepareImportanceSampling(self.istore)
        self.sampler.PrepareScoring(self.scorer)
        self.sampler.Configure()






    def run(self, nruns):
        if not "sampler" in dir(self):
            print "call \"setupSampler\" first"
            return
        # this function class BeamOn(self.nevents) nruns times
        # and fills a histogram with the number of tracks
        # entered the 4th cell in the self.CellScorerList ("p_D2")
        # at the end it prints a table of the scores

        
        # get run manager and run 100 events
        
        r = self.base.GetRunManager()
        tinold=self.CellScorerList[17].GetCellScoreValues().fSumTracksEntering
        # run nruns runs with self.nevents events
        # each and fill a histogram after each run
        for i in range(nruns):
            r.BeamOn(self.nevents)
            tin = self.CellScorerList[17].GetCellScoreValues().fSumTracksEntering
            v = tin-tinold
            self.messurements.append(v)
            tinold = tin
            print "run number:", i+1, "messured: ", v

        # create a G4ScoreTable
        st = G4ScoreTable(self.istore)
        # get the result as table in a string
        s = st.Write(self.cs_store.GetMapGeometryCellCellScorer())
        print s


