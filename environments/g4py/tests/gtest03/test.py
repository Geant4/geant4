#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   gtest02
#   - test for using site-module packages
# ==================================================================
from Geant4 import *
import NISTmaterials
from  EZsim import EZgeom
from EZsim.EZgeom import G4EzVolume
import EMSTDpl
import ParticleGun

# ==================================================================
# intialize
# ==================================================================
def Configure():
  # ------------------------------------------------------------------
  # setup for materials
  # ------------------------------------------------------------------
  # simple materials for Qgeom
  NISTmaterials.Construct()

  # ------------------------------------------------------------------
  # setup for geometry
  # ------------------------------------------------------------------
  #Qgeom.Construct()
  EZgeom.Construct()  # initialize

  # ------------------------------------------------------------------
  # setup for physics list
  # ------------------------------------------------------------------
  EMSTDpl.Construct()

  # ------------------------------------------------------------------
  # setup for primary generator action
  # ------------------------------------------------------------------
  ParticleGun.Construct()
  gControlExecute("gun.mac")

# ==================================================================
# constructing geometry
# ==================================================================
def ConstructGeom():
  print "* Constructing geometry..."
  # reset world material
  air= G4Material.GetMaterial("G4_AIR")
  EZgeom.SetWorldMaterial(air)

  # target
  global target
  target= G4EzVolume("Target")
  au= G4Material.GetMaterial("G4_Au")
  target.CreateTubeVolume(au, 0., 1.*cm, 1.*mm)
  target.PlaceIt(G4ThreeVector(0.,0.,-10.*cm))

  # dummy box
  global detector_box, detector_box_pv
  detector_box= G4EzVolume("DetectorBox")
  detector_box.CreateBoxVolume(air, 20.*cm, 20.*cm, 40.*cm)
  detector_box_pv= detector_box.PlaceIt(G4ThreeVector(0.,0.,20.*cm))

  # calorimeter
  global cal
  cal= G4EzVolume("Calorimeter")
  nai= G4Material.GetMaterial("G4_SODIUM_IODIDE")
  cal.CreateBoxVolume(nai, 5.*cm, 5.*cm, 30.*cm)
  dd= 5.*cm
  for ical in range(-1, 2):
    calPos= G4ThreeVector(dd*ical, 0., 0.)
    print calPos
    cal.PlaceIt(calPos, ical+1, detector_box)


# ==================================================================
# main
# ==================================================================
# ------------------------------------------------------------------
# randum number
# ------------------------------------------------------------------
rand_engine= Ranlux64Engine()
HepRandom.setTheEngine(rand_engine)
HepRandom.setTheSeed(20050830L)

# setup...
Configure()
ConstructGeom()

# ------------------------------------------------------------------
# go...
# ------------------------------------------------------------------
gRunManager.Initialize()

# visualization
gControlExecute("vis.mac")

# beamOn
#gRunManager.BeamOn(3)

