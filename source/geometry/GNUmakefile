# $Id: GNUmakefile,v 1.5 2003-05-28 07:54:30 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for geometry.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------

name := G4geometry

SUBDIRS  = management volumes magneticfield verification biasing
SUBDIRS += solids/CSG solids/specific solids/Boolean solids/BREPS
SUBLIBS  = G4geometrymng G4volumes G4magneticfield G4geomver G4geombias
SUBLIBS += G4csg G4specsolids G4geomBoolean G4brep
ifdef G4USE_STEP
  SUBDIRS += solids/STEP solids/STEPinterface
  SUBLIBS += G4step G4stepinterface
endif

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/globlib.gmk
