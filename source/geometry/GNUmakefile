# $Id: GNUmakefile,v 1.1 1999-01-07 16:06:59 gunter Exp $
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
