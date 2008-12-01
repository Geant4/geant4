#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   gtest02
#   - test for using site-module packages
# ==================================================================
from Geant4 import *
import g4py.NISTmaterials
import g4py.ezgeom
from g4py.ezgeom import G4EzVolume
import g4py.EMSTDpl
import g4py.ParticleGun

# ==================================================================
# intialize
# ==================================================================
def Configure():
  # ------------------------------------------------------------------
  # setup for materials
  # ------------------------------------------------------------------
  # simple materials for Qgeom
  g4py.NISTmaterials.Construct()

  # ------------------------------------------------------------------
  # setup for geometry
  # ------------------------------------------------------------------
  #g4py.Qgeom.Construct()
  g4py.ezgeom.Construct()  # initialize

  # ------------------------------------------------------------------
  # setup for physics list
  # ------------------------------------------------------------------
  g4py.EMSTDpl.Construct()

  # ------------------------------------------------------------------
  # setup for primary generator action
  # ------------------------------------------------------------------
  g4py.ParticleGun.Construct()
  gControlExecute("gun.mac")


# ==================================================================
# constructing geometry
# ==================================================================
def ConstructGeom():
  print "* Constructing geometry..."
  # reset world material
  air= G4Material.GetMaterial("G4_AIR")
  g4py.ezgeom.SetWorldMaterial(air)

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

