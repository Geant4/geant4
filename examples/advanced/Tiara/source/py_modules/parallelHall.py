# $Id: parallelHall.py,v 1.4 2003/06/20 12:41:06 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-09-00 $
# -------------------------------------------------------------------
#
import G4Kernel
import CLHEP
import string

class ParallelHall(object):
    def __init__(self, tiaraSpecs):
        self.halfWidth = tiaraSpecs.dimensions.worldHalfWidth + \
                     10*CLHEP.cm
        self.halfLength = tiaraSpecs.dimensions.worldHalfLength + \
                      10*CLHEP.cm
        self.hallSolid = G4Kernel.G4Box("parallelBox", 
                                        self.halfWidth,
                                        self.halfWidth,
                                        self.halfLength)
        vacuum = tiaraSpecs.materials.GetMaterial("vacuum")
        self.logHall = G4Kernel.G4LogicalVolume(self.hallSolid,
                                                vacuum,
                                                "paralleleLog")
        rot = G4Kernel.G4RotationMatrix()
        self.worldVolume = G4Kernel.\
                           G4PVPlacement(rot,
                                         CLHEP.Hep3Vector(0, 0, 0),
                                         self.logHall,
                                         "ParallelHall");
        self.geoCells = []

    def getWorldVolume(self):
        return self.worldVolume

    def placeCells(self, arrPosLogVol):
        for ele in arrPosLogVol:
            phys = self.placeOneCell(ele)
            self.geoCells.append(G4Kernel.G4GeometryCell(phys, 0))

    def placeOneCell(self, pLog, name = ""):
        rot = G4Kernel.G4RotationMatrix()
        z = int(pLog.pos)
        zstr = "%(z)d" % vars()
        zstr = string.rjust(zstr, 5)
        zstr = string.replace(zstr,' ','0')

        if not name:
            name = "cell_z_" + zstr
        vPhys = G4Kernel.\
                G4PVPlacement(rot,
                              CLHEP.Hep3Vector(0, 0, pLog.pos),
                              pLog.log,
                              name,
                              self.logHall);
        return vPhys

    def getGeometryCells(self):
        return self.geoCells


