# $Id: GNUmakefile,v 1.3 2003/09/21 19:38:49 kurasige Exp $
# --------------------------------------------------------------
# GNUmakefile for particles library.  Gabriele Cosmo, 4/4/97.
# --------------------------------------------------------------

name := G4particles

SUBDIRS = management bosons leptons shortlived
SUBDIRS += hadrons/mesons hadrons/barions hadrons/ions
SUBDIRS += utils
SUBLIBS = G4partman G4bosons G4leptons G4shortlived
SUBLIBS += G4mesons G4baryons G4ions
SUBLIBS += G4partutils
ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
