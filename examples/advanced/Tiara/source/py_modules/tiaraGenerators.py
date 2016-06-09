# $Id: tiaraGenerators.py,v 1.2 2003/06/16 17:06:44 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-06-00 $
# -------------------------------------------------------------------
#
import string
import tiaraSpecifications
import Tiara
import os

tiara_dir = os.environ["TIARA_BASE"]

class TiaraPrimaryGenerator(object):
    def __init__(self, tiaraSpecs):
        self.name = "preIntegratedSource"
        self.eSamp = Tiara.TiaraSampledEnergy(\
            tiaraSpecs.experiment.energy,
            tiaraSpecs.experiment.minNeutronEnergyCut,
            tiara_dir + "/data/expDataConverted/source.xml",
            "_v3")
        self.tally = Tiara.TiaraTally()
        self.tally.setBinEdges(tiaraSpecifications.\
            sourceTallyEdges[tiaraSpecs.experiment.energy])
        
        self.directionGenerator = Tiara.TiaraIsotropicDirections(\
            0,
            tiaraSpecs.dimensions)
        self.primGen = Tiara.TiaraPrimaryGeneratorAction(\
            self.eSamp,
            self.directionGenerator,
            self.tally,
            tiaraSpecs.dimensions)
        return


class TiaraDPSEnergyGenerator(object):
    def __init__(self, tiaraSpecs, xmlName):
        self.Name = "dpsSource"

        self.eSamp = Tiara.\
                     TiaraDPSSampledEnergy(tiaraSpecs.
                                           experiment.energy,
                                           tiaraSpecs.
                                           experiment.minNeutronEnergyCut,
                                           xmlName,
                                           "")
        

        self.tally = Tiara.TiaraTally()

        if string.find(xmlName,"other") > -1:
            if tiaraSpecs.experiment.energy == "43":
                self.tally.setBinEdges([0,37.5, 43,50])
            else:
                if tiaraSpecs.experiment.energy == "68":
                    self.tally.setBinEdges([0,61, 69,80])
        else:
            self.tally.setBinEdges(tiaraSpecifications.\
                                   sourceTallyEdges[tiaraSpecs.
                                                    experiment.energy])

        self.directionGenerator = Tiara.TiaraIsotropicDirections(\
            0,
            tiaraSpecs.dimensions)
        self.primGen = Tiara.TiaraPrimaryGeneratorAction(\
            self.eSamp,
            self.directionGenerator,
            self.tally,
            tiaraSpecs.dimensions)
        return
        
class FixedEnergyPrimaryGenerator(object):
    def __init__(self, tiaraSpecs):
        self.name = "fixedSource"
        self.eSamp = Tiara.TiaraFixedEnergyGenerator(\
            float(tiaraSpecs.experiment.energy))
        self.tally = Tiara.TiaraTally()
        self.tally.setBinEdges(tiaraSpecifications.\
                               sourceTallyEdges[tiaraSpecs.
                                                experiment.energy])

        self.directionGenerator = Tiara.TiaraIsotropicDirections(\
            0,
            tiaraSpecs.dimensions)
        self.primGen = Tiara.TiaraPrimaryGeneratorAction(\
            self.eSamp,
            self.directionGenerator,
            self.tally,
            tiaraSpecs.dimensions)
        return



