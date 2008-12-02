#!/usr/bin/python
# ==================================================================
# An example of writing a GDML file
#
# ==================================================================
from Geant4 import *
import g4py.Qmaterials, g4py.Qgeom
import g4py.ExN01pl, g4py.ParticleGun

# ==================================================================
# main
# ==================================================================
# set geometry
g4py.Qmaterials.Construct()
g4py.Qgeom.Construct()

# minimal physics list
g4py.ExN01pl.Construct()

# set primary generator action
g4py.ParticleGun.Construct()

# initialize
gRunManager.Initialize()

# visualization
gApplyUICommand("/vis/open OGLIX")
gApplyUICommand("/vis/scene/create")
gApplyUICommand("/vis/scene/add/volume")
gApplyUICommand("/vis/sceneHandler/attach")
gApplyUICommand("/vis/viewer/set/viewpointThetaPhi 90. -90.")

# write to a GDML file
print "\n*** write to a GDML file..."
navigator= gTransportationManager.GetNavigatorForTracking()
world_volume= navigator.GetWorldVolume()

gdml_parser = G4GDMLParser()
gdml_parser.Write("qgeom.gdml", world_volume)


