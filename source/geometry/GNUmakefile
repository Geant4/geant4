# $Id: GNUmakefile,v 1.4 2002-03-25 17:09:16 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for geometry.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------

name := G4geometry

SUBDIRS  = management volumes magneticfield verification biasing
SUBDIRS += solids/CSG solids/specific solids/Boolean solids/BREPS solids/STEP solids/STEPinterface 
SUBLIBS  = G4geometrymng G4volumes G4magneticfield G4geomver G4geombias
SUBLIBS += G4csg G4specsolids G4geomBoolean G4brep G4step G4stepinterface

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/globlib.gmk
