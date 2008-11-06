# $Id: GNUmakefile,v 1.5 2008-11-06 13:17:36 kurasige Exp $
# --------------------------------------------------------------
# GNUmakefile for particles library.  Gabriele Cosmo, 4/4/97.
# --------------------------------------------------------------

name := G4particles

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS  = libG4geometry.lib libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

SUBDIRS = management bosons leptons shortlived
SUBDIRS += hadrons/mesons hadrons/barions hadrons/ions
SUBDIRS += adjoint utils 
SUBLIBS = G4partman G4bosons G4leptons G4shortlived
SUBLIBS += G4mesons G4baryons G4ions
SUBLIBS += G4partadj G4partutils

include $(G4INSTALL)/config/globlib.gmk
