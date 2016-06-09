# $Id: GNUmakefile,v 1.11 2003/10/08 14:39:19 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for geometry.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------

name := G4geometry

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

SUBDIRS  = management volumes magneticfield navigation biasing divisions
SUBDIRS += solids/CSG solids/specific solids/Boolean solids/BREPS
SUBLIBS  = G4geometrymng G4volumes G4magneticfield
SUBLIBS += G4navigation G4geombias G4geomdivision
SUBLIBS += G4csg G4specsolids G4geomBoolean G4brep

include $(G4INSTALL)/config/globlib.gmk
