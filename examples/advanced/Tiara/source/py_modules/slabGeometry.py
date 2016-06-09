# $Id: slabGeometry.py,v 1.3 2003/06/20 12:41:07 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-09-01 $
# -------------------------------------------------------------------
#
import math
import G4Kernel
import parallelHall
import posLog

class SlabedGeometry(object):
    def __init__(self, tiaraSpecs, cellWidth, parallelGeo):
        self.arrPosLogVol = []
        self.normalCellSolid = None
        self.logNormCell = None
        self.lastCellSolid = None
        self.logLastCell = None
        self.createCells(tiaraSpecs, cellWidth, parallelGeo)

    def createCells(self, tiaraSpecs, cellWidth, parallelGeo):
        print "SlabedGeometry::createCells:"
        print "halfLength", parallelGeo.halfLength,\
              "halfWidth", parallelGeo.halfWidth
        nCells = int(tiaraSpecs.experiment.shieldWidth / cellWidth)
        print "cellWidth", cellWidth
        print "nCells", nCells
        zCellStart = tiaraSpecs.dimensions.targetPosZ + \
                     tiaraSpecs.dimensions.distTargetExperiment + \
                     tiaraSpecs.experiment.colWidth
        print "zCellStart", zCellStart
        zCellRegion = parallelGeo.halfLength - zCellStart
        print "zCellRegion", zCellRegion
        lengthLastCell = zCellRegion - nCells * cellWidth
        print "lengthLastCell", lengthLastCell
        print "sum cell length:", lengthLastCell + nCells * cellWidth
        vacuum = tiaraSpecs.materials.GetMaterial("vacuum")
        self.normalCellSolid = G4Kernel.G4Box("cellBox", 
                                              parallelGeo.halfWidth,
                                              parallelGeo.halfWidth,
                                              cellWidth/2)
        print "normal box:",parallelGeo.halfWidth,cellWidth/2
        self.logNormCell = G4Kernel.G4LogicalVolume(self.normalCellSolid,
                                                    vacuum,
                                                    "cellLog")

        zPos = zCellStart - 0.5 * cellWidth
        for i in range(nCells):
            zPos += cellWidth
            self.arrPosLogVol.append(posLog.PosLog(zPos, self.logNormCell))
            print "zPos", zPos
            
            
        self.lastCellSolid = G4Kernel.G4Box("lastCellBox", 
                                            parallelGeo.halfWidth,
                                            parallelGeo.halfWidth,
                                            lengthLastCell/2)
        self.logLastCell = G4Kernel.G4LogicalVolume(self.lastCellSolid,
                                                    vacuum,
                                                    "lastCellLog")

        print "last cell:", parallelGeo.halfWidth, lengthLastCell/2

        zPos += cellWidth/2+lengthLastCell/2
        self.arrPosLogVol.append(posLog.PosLog(zPos, self.logLastCell))
        print "last cell zPos", zPos
            


        
    def getArrPosLogVol(self):
        return self.arrPosLogVol



class SlabedImportanceGeometry(object):
    def __init__(self, tiaraSpecs, cellWidth, impBase, parallelGeo = None):
        self.parallelGeo = parallelGeo
        self.tiaraSpecs = tiaraSpecs
        self.cellWidth = cellWidth
        self.iStore = None
        self.geometryCells = []
        self.base = impBase
        cellwidth_cm = cellWidth / CLHEP.cm
        self.nameExt = "-cellWidth_%(cellwidth_cm)d" %vars()
        self.buildParallelGeometry()
        self.setImportances()
        
    def buildParallelGeometry(self):
        self.parallelGeo = parallelHall.ParallelHall(self.tiaraSpecs)

        self.slabedGeo = SlabedGeometry(self.tiaraSpecs,
                                        self.cellWidth,
                                        self.parallelGeo)
        self.parallelGeo.placeCells(self.slabedGeo.getArrPosLogVol())
        self.iStore = G4Kernel.G4IStore(self.parallelGeo.getWorldVolume())
        self.geometryCells = self.parallelGeo.getGeometryCells()
        
    def setImportances(self):
        worldCell = G4Kernel.G4GeometryCell(self.parallelGeo.\
                                            getWorldVolume(), 0)
        self.iStore.AddImportanceGeometryCell(1, worldCell)
        nCells = len(self.geometryCells)
        for i in range(nCells-1):
            cell = self.geometryCells[i]
            importance = math.pow(self.base, i)
            print "i=", importance
            self.iStore.AddImportanceGeometryCell(importance, cell)
        lastCell = self.geometryCells[nCells-1]
        importance = math.pow(self.base, nCells-2)
        print "last cells i=", importance
        self.iStore.AddImportanceGeometryCell(importance, lastCell)
        

            
    def getWorldVolume(self):
        return self.parallelGeo.getWorldVolume()
    
    def getImportanceStore(self):
        return self.iStore
