# $Id: tiaraSpecifications.py,v 1.3 2003/06/20 12:41:07 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-09-00 $
# -------------------------------------------------------------------
#
import CLHEP
import G4Kernel

ColWidth = {"43":40*CLHEP.cm,
            "68":80*CLHEP.cm}

BinEdgesScinti = {"43":range(4,45),
                  "68":range(6,45) + range(46,71,2)}

BinEdgesBonner = {"43":       [4.500E+07,
                               3.500E+07,
                               2.750E+07,
                               2.250E+07,
                               1.750E+07,
                               1.350E+07,
                               1.000E+07,
                               6.700E+06,
                               4.490E+06,
                               3.010E+06,
                               2.020E+06,
                               1.350E+06,
                               9.070E+05,
                               4.980E+05,
                               2.240E+05,
                               8.650E+04,
                               1.500E+04,
                               3.350E+03,
                               4.540E+02,
                               2.260E+01,
                               5.040E+00,
                               1.120E+00,
                               4.140E-01,
                               1.000E-04],
                  "68":       [ 8.000E+07,
                                6.500E+07,
                                5.500E+07,
                                4.500E+07,
                                3.500E+07,
                                2.750E+07,
                                2.250E+07,
                                1.750E+07,
                                1.350E+07,
                                1.000E+07,
                                6.700E+06,
                                4.490E+06,
                                3.010E+06,
                                2.020E+06,
                                1.350E+06,
                                9.070E+05,
                                4.980E+05,
                                2.240E+05,
                                8.650E+04,
                                1.500E+04,
                                3.350E+03,
                                4.540E+02,
                                2.260E+01,
                                5.040E+00,
                                1.120E+00,
                                4.140E-01,
                                1.000E-04]}


tallyBinEdges = {}
tallyBinEdges["43"] = [0,10,35,45]
tallyBinEdges["68"] = [0,10,60,70]

sourceTallyEdges = {}
sourceTallyEdges["43"] = [0,36.3, 45.5,50]
sourceTallyEdges["68"] = [0,60.8, 72.5,80]

class Experiment(object):
    def __init__(self,
                 energy,
                 minNeutronEnergyCut,
                 particleCut,
                 shieldWidth,
                 shieldMaterial):
        self.energy = "%(energy)d" % vars()
        self.minNeutronEnergyCut = minNeutronEnergyCut
        self.particleCut = particleCut
        self.shieldWidth = shieldWidth
        self.shieldMaterial = shieldMaterial
        self.binEdgesBonner = BinEdgesBonner[self.energy]
        self.binEdgesScinti = BinEdgesScinti[self.energy]
        self.colWidth = 0
        if self.shieldWidth <= 50.0*CLHEP.cm:
            self.colWidth = ColWidth[self.energy]


class Specifications(object):
    def __init__(self,dimensions, experiment, materials):
        self.dimensions = dimensions
        self.experiment = experiment
        self.materials = materials

