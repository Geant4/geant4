# $Id: GNUmakefile,v 1.2 2000-04-07 13:07:00 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for geometry.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------

name := G4geometry

SUBDIRS  = management volumes magneticfield
SUBDIRS += solids/CSG solids/specific solids/Boolean solids/BREPS solids/STEP solids/STEPinterface 
SUBLIBS  = G4geometrymng G4volumes G4magneticfield
SUBLIBS += G4csg G4specsolids G4geomBoolean G4brep G4step G4stepinterface

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/globlib.gmk
