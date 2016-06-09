# $Id: variableGeometry.py,v 1.2 2003/06/16 17:06:45 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-09-01 $
# -------------------------------------------------------------------
#
import G4Kernel
import parallelHall
import posLog

class VariableSlobedGeometry(object):
    def __init__(self,\
                 tiaraSpecs,
                 cellSizeImportanceList,
                 parallelGeo):
        self.arrPosLogVol = []
        self.normalCellSolids = []
        self.logNormCells = []
        self.lastCellSolid = None
        self.logLastCell = None
        self.createCells(tiaraSpecs,
                         cellSizeImportanceList, parallelGeo)

    def createCells(self, tiaraSpecs,
                    cellSizeImportanceList, parallelGeo):
        nImps = len(cellSizeImportanceList)
        zCellStart = tiaraSpecs.dimensions.targetPosZ + \
                     tiaraSpecs.dimensions.distTargetExperiment + \
                     tiaraSpecs.experiment.colWidth
        zCellRegion = parallelGeo.halfLength - zCellStart
        print "zCellStart", zCellStart, "zCellRegion", zCellRegion
        vacuum = tiaraSpecs.materials.GetMaterial("vacuum")
        for i in range(nImps):
            ncsName = "cellBox_%(i)d" % vars()
            print "i", i, "cswidth", cellSizeImportanceList[i]["width"]
            ncs = G4Kernel.G4Box(ncsName,
                                parallelGeo.halfWidth,
                                parallelGeo.halfWidth,
                                0.5 * cellSizeImportanceList[i]["width"])
            self.normalCellSolids.append(ncs)
            nclName = "cellLog_%(i)d" % vars()
            ncl = G4Kernel.G4LogicalVolume(ncs,
                                           vacuum,
                                           nclName)
            self.logNormCells.append(ncl)

        zPos = zCellStart
        for i in range(nImps):
            zPos += 0.5 * cellSizeImportanceList[i]["width"]
            print "zPos", zPos, "cellNum", i
            self.arrPosLogVol.append(posLog.PosLog(zPos, self.logNormCells[i]))
            zPos += 0.5 * cellSizeImportanceList[i]["width"]

            
        allCellWidth = 0.0
        for i in range(nImps):
            allCellWidth+=cellSizeImportanceList[i]["width"]
        print "allCellWidth", allCellWidth
        lengthLastCell = zCellRegion - allCellWidth
            
        self.lastCellSolid = G4Kernel.G4Box("lastCellBox", 
                                            parallelGeo.halfWidth,
                                            parallelGeo.halfWidth,
                                            lengthLastCell/2)
        self.logLastCell = G4Kernel.G4LogicalVolume(self.lastCellSolid,
                                                    vacuum,
                                                    "lastCellLog")

        zPos += lengthLastCell/2
        print "last cell zPos", zPos, "length", lengthLastCell
        self.arrPosLogVol.append(posLog.PosLog(zPos, self.logLastCell))

        
    def getArrPosLogVol(self):
        return self.arrPosLogVol
        

class VariableImpSlabGeometry(object):
    def __init__(self, tiaraSpecs,
                 parallelGeo = None):
        self.tiaraSpecs = tiaraSpecs
        self.cellSizeImportanceList = []
        self.parallelGeo = parallelGeo
        self.iStore = None
        self.geometryCells = []
        self.base = 0.0
        self.nameExt = "-variableCellWidth"


    def addCellImportance(self, width, faktor):
        self.cellSizeImportanceList.append(
            {"width":width, "faktor":faktor})

        
    def construct(self):
        nImps = len(self.cellSizeImportanceList)
        if nImps>0:
            for i in range(nImps):
                self.base+=1.0*self.cellSizeImportanceList[i]["faktor"]
            self.base/=nImps
        self.createParallelGeometry(self.tiaraSpecs)
        self.setImportances()
        

    def createParallelGeometry(self, tiaraSpecs):
        self.parallelGeo = parallelHall.ParallelHall(tiaraSpecs)
        self.slobedGeo = VariableSlobedGeometry(\
            tiaraSpecs,
            self.cellSizeImportanceList,
            self.parallelGeo)
        self.parallelGeo.placeCells(self.slobedGeo.getArrPosLogVol())
        self.iStore = G4Kernel.G4IStore(self.parallelGeo.getWorldVolume())
        self.geometryCells = self.parallelGeo.getGeometryCells()
        
    def setImportances(self):
        worldCell = G4Kernel.G4GeometryCell(self.parallelGeo.\
                                            getWorldVolume(), 0)
        self.iStore.AddImportanceGeometryCell(1, worldCell)
        nCells = len(self.geometryCells)
        nImps = len(self.cellSizeImportanceList)
        if ((nCells - 1)  != nImps):
            print "VariableImpSlabGeometry: ERROR importances and cells don't mach!", nCells, nImps
        importance = 1
        for i in range(nImps):
            cell = self.geometryCells[i]
            importance*=self.cellSizeImportanceList[i]["faktor"]
            print "i=", importance
            self.iStore.AddImportanceGeometryCell(importance, cell)

        lastCell = self.geometryCells[nCells-1]
        print "last cells i=", importance
        self.iStore.AddImportanceGeometryCell(importance, lastCell)

    def getWorldVolume(self):
        return self.parallelGeo.getWorldVolume()

    def getImportanceStore(self):
        return self.iStore
