# $Id: GNUmakefile,v 1.3 2001-10-17 12:59:48 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for geometry.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------

name := G4geometry

SUBDIRS  = management volumes magneticfield verification
SUBDIRS += solids/CSG solids/specific solids/Boolean solids/BREPS solids/STEP solids/STEPinterface 
SUBLIBS  = G4geometrymng G4volumes G4magneticfield G4geomver
SUBLIBS += G4csg G4specsolids G4geomBoolean G4brep G4step G4stepinterface

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/globlib.gmk
