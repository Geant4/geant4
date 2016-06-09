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

  # phantom
  global phantom
  phantom= G4EzVolume("DetectorBox")
  water= G4Material.GetMaterial("G4_WATER")
  phantom.CreateBoxVolume(water, 40.*cm, 40.*cm, 50.*cm)
  phantom.PlaceIt(G4ThreeVector(0.,0.,20.*cm))
  vsize=phantom.VoxelizeIt(100, 100, 100)
  print "voxel size=", vsize

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

