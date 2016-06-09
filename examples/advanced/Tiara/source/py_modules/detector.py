# $Id: detector.py,v 1.2 2003/06/16 17:06:44 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-05-02 $
# -------------------------------------------------------------------
#

# A detector is either the of the "simple" or "ring" type.
# the volumina of the "simple" detectors are the same at all
# positions (00, 20, 40 cm off beam axis). The "ring" type detectors
# have different volumina at the different positions.
# For convinience a class Detector is provided to link a detector
# position and it's volume.

d00Volume = 1608.8   # volume of the 12.9 cm x 12.9 cm cyinder
d20Volume = 20268.3  # volume of the ring at 20 cm
d40Volume = 40536.6  # volume of the ring at 40 cm

# map of distances to volumina for the 12.9 cm x 12.9 cm detectors
SimpleDetectorVolume = {"00":d00Volume, "20":d00Volume, "40":d00Volume}
# map of distances to volumina for ring detectors
RingDetectorVolume = {"00":d00Volume, "20":d20Volume, "40":d40Volume}

# map of detector type to distance-volumina map
DetectorVolume = {"simple":SimpleDetectorVolume,
                  "ring"  :RingDetectorVolume}

class Detector(object):
    "Distance and volume of a Detector."
    def __init__(self, dist, detType):
        self.name = dist
        self.volume = DetectorVolume[detType][dist]





# function to calculate the scale factor

def detScale(ngen, energy, coli):
    """Detrmine scale.

    Determine scale according to the number of generated neutrons,
    the proton beam energy. If coli == 1, the flux at the
    colimator exit is used, else the flux at 401 cm is used.
    """
    Asrc = 93.31
    Fexp = 0
    if energy == "43":
        if coli == 1:
            Fexp = 1.76E+04
        else:
            Fexp = 1.94E+04
    else:
        if energy == "68":
            if coli == 1:
                Fexp = 2.04E+04
            else:
                Fexp = 2.46E+04

    S = Fexp * Asrc / ngen

    return S


