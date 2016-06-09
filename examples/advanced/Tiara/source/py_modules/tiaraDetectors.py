# $Id: tiaraDetectors.py,v 1.4 2004/06/09 15:04:36 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-08-00 $
# -------------------------------------------------------------------
#
import CLHEP
import G4Kernel

class ScoreDetector(object):
    def __init__(self, name):
        self.name = name
        self.phys = ""
        self.scorer = ""

class DetectorSlab(object):
    def __init__(self):
        pass
    def createScoreDetectors(self, tiaraHall):
        scoreDets = []
        scoreDets.append(ScoreDetector("detectorSlab"))
        scoreDets[0].phys = tiaraHall.AddDetectorSlab(scoreDets[0].name)
        print "+++ DetectorSlab: created slab detector."
        return scoreDets

class ThreeZylindricDetectors(object):
    def __init__(self):
        pass
    def createScoreDetectors(self, tiaraHall):
        scoreDets = []
        scoreDets.append(ScoreDetector("detector_00"))
        scoreDets.append(ScoreDetector("detector_20"))
        scoreDets.append(ScoreDetector("detector_40"))
        dist = 0.0
        for det in scoreDets:
            if dist > 0.0:
                det.phys = tiaraHall.AddPhysicalRingDetector(dist,
                                                             det.name)
            else:
                det.phys = tiaraHall.AddPhysicalDetector(dist, det.name)
            dist += 20*CLHEP.cm
        print "+++ ThreeZylindricDetectors: created 3 zylindric detectors"
        return scoreDets
