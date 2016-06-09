#!/usr/bin/python3
# ==================================================================
# python script for Geant4Py test
#
#   gtest02
#   - test for using site-module packages
# ==================================================================
from Geant4 import *
import g4py.Qmaterials, g4py.NISTmaterials
import g4py.Qgeom, g4py.ExN01geom, g4py.ExN03geom
import g4py.ExN01pl, g4py.EMSTDpl
import g4py.ParticleGun, g4py.MedicalBeam

# ==================================================================
# user setup
# ==================================================================

# ------------------------------------------------------------------
# Setup-0 (Q)
# ------------------------------------------------------------------
def Setup0():
  # simple materials for Qgeom
  g4py.Qmaterials.Construct()

  # NIST materials
  #g4py.NISTmaterials.Construct()
  
  # normal way for constructing user geometry
  #qDC= g4py.Qgeom.QDetectorConstruction()
  #gRunManager.SetUserInitialization(qDC)

  # 2nd way, short-cut way
  g4py.Qgeom.Construct()

  # primary
  global primary_position, primary_direction
  primary_position= G4ThreeVector(0.,0., -14.9*cm)
  primary_direction= G4ThreeVector(0.2, 0., 1.)


# ------------------------------------------------------------------
# Setup-1 (ExampleN01)
# ------------------------------------------------------------------
def Setup1():
  g4py.ExN01geom.Construct()

  global primary_position, primary_direction
  primary_position= G4ThreeVector(-2.5*m, 0., 0.)
  primary_direction= G4ThreeVector(1., 0., 0.)


# ------------------------------------------------------------------
# Setup-3 (ExampleN03)
# ------------------------------------------------------------------
def Setup3():
  #exN03geom= g4py.ExN03geom.ExN03DetectorConstruction()
  #gRunManager.SetUserInitialization(exN03geom)

  g4py.ExN03geom.Construct()

  global primary_position, primary_direction
  primary_position= G4ThreeVector(-1.*m, 0., 0.)
  primary_direction= G4ThreeVector(1., 0., 0.)


# ==================================================================
# main
# ==================================================================
# ------------------------------------------------------------------
# randum number
# ------------------------------------------------------------------
rand_engine= Ranlux64Engine()
HepRandom.setTheEngine(rand_engine)
HepRandom.setTheSeed(20050830)

# ------------------------------------------------------------------
# user setup
# ------------------------------------------------------------------
Setup0()
#Setup1()
#Setup3()


# ------------------------------------------------------------------
# setup for physics list
# ------------------------------------------------------------------
# normal way for constructing user physics list
#exN01PL= ExN01PhysicsList.ExN01PhysicsList()
#gRunManager.SetUserInitialization(exN01PL)

# 2nd way, short-cut way
# geantino + transportation
#g4py.ExN01pl.Construct()

# electron/gamma standard EM
g4py.EMSTDpl.Construct()

# ------------------------------------------------------------------
# setup for primary generator action
# ------------------------------------------------------------------
# ------------
# Particle Gun
# ------------
# normal way for constructing user physics list
#pgPGA= g4py.ParticleGun.ParticleGunAction()
#gRunManager.SetUserAction(pgPGA)
#pg= pgPGA.GetParticleGun()

# 2nd way, short-cut way
pg= g4py.ParticleGun.Construct()

# set parameters of particle gun
pg.SetParticleByName("e-")
pg.SetParticleEnergy(300.*MeV)
pg.SetParticlePosition(primary_position)
pg.SetParticleMomentumDirection(primary_direction)

# ------------
# Medical Beam
# ------------
#beam= g4py.MedicalBeam.Construct()

# ------------------------------------------------------------------
# go...
# ------------------------------------------------------------------
gRunManager.Initialize()

# visualization
gApplyUICommand("/control/execute vis.mac")

# beamOn
#gRunManager.BeamOn(3)

