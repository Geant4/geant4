# $Id: GNUmakefile,v 2.2 1998/10/30 18:55:38 japost Exp $
# --------------------------------------------------------------
# GNUmakefile for geometry.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------

name := G4geometry

SUBDIRS = management volumes solids/CSG solids/Boolean solids/BREPS solids/STEP solids/STEPinterface magneticfield 
SUBLIBS = G4geometrymng G4volumes G4csg G4geomBoolean G4brep G4step G4stepinterface G4magneticfield

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/globlib.gmk
