# $Id: GNUmakefile,v 1.1 1999-01-07 16:09:46 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for particles library.  Gabriele Cosmo, 4/4/97.
# --------------------------------------------------------------

name := G4particles

SUBDIRS = management bosons leptons shortlived
SUBDIRS += hadrons/mesons hadrons/barions hadrons/ions
SUBLIBS = G4partman G4bosons G4leptons G4shortlived
SUBLIBS += G4mesons G4barions G4ions

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
