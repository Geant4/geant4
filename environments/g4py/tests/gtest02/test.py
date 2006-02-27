#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   gtest02
#   - test for using site-module packages
# ==================================================================
from Geant4 import *
import Qmaterials, NISTmaterials
import Qgeom, ExN01geom, ExN03geom
import ExN01pl, EMSTDpl
import ParticleGun, MedicalBeam

# ==================================================================
# main
# ==================================================================
# ------------------------------------------------------------------
# randum number
# ------------------------------------------------------------------
rand_engine= Ranlux64Engine()
HepRandom.setTheEngine(rand_engine)
HepRandom.setTheSeed(20050830L)

# ------------------------------------------------------------------
# setup for materials
# ------------------------------------------------------------------
# simple materials for Qgeom
#Qmaterials.Construct()

# NIST materials
#NISTmaterials.Construct()

# ------------------------------------------------------------------
# setup for geometry
# ------------------------------------------------------------------
# normal way for constructing user geometry
#qDC= Qgeom.QDetectorConstruction()
#gRunManager.SetUserInitialization(qDC)

exN03geom= ExN03geom.ExN03DetectorConstruction()
gRunManager.SetUserInitialization(exN03geom)

# 2nd way, short-cut way
#Qgeom.Construct()
#ExN01geom.Construct()
#ExN03geom.Construct()

# ------------------------------------------------------------------
# setup for physics list
# ------------------------------------------------------------------
# normal way for constructing user physics list
#exN01PL= ExN01pl.ExN01PhysicsList()
#gRunManager.SetUserInitialization(exN01PL)

# 2nd way, short-cut way
#ExN01pl.Construct()
EMSTDpl.Construct()

# ------------------------------------------------------------------
# setup for primary generator action
# ------------------------------------------------------------------
# normal way for constructing user physics list
#pgPGA= ParticleGun.ParticleGunAction()
#gRunManager.SetUserAction(pgPGA)
#pg= pgPGA.GetParticleGun()

# 2nd way, short-cut way
pg= ParticleGun.Construct()

# set parameters of particle gun
pg.SetParticleByName("e-")
pg.SetParticleEnergy(200.*MeV)
pg.SetParticlePosition(G4ThreeVector(0.,0.,-14.9)*cm)
pg.SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.))

# medical beam
#beam= MedicalBeam.Construct()

# ------------------------------------------------------------------
# go...
# ------------------------------------------------------------------
gRunManager.Initialize()

# visualization
gApplyUICommand("/control/execute vis.mac")

# beamOn
#gRunManager.BeamOn(3)

