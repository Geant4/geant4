#!/usr/bin/python
# ==================================================================
# An example of writing a GDML file
#
# ==================================================================
from Geant4 import *
from Geant4.G4gdml import *

import Qmaterials, Qgeom
import ExN01pl, ParticleGun

# ==================================================================
# main
# ==================================================================
# set geometry
Qmaterials.Construct()
Qgeom.Construct()

# minimal physics list
ExN01pl.Construct()

# set primary generator action
ParticleGun.Construct()

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
gdml_writer= G4GDMLWriter("GDMLSchema/gdml.xsd", "qgeom.gdml")
gdml_writer.DumpGeometryInfo(world_volume)

