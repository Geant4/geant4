# $Id$
# ----------------------------------------------------------------
# GNUmakefile for hadronic processes library.  G.Folger 10-Dec-97.
# ----------------------------------------------------------------

name := G4biasing

SUBDIRS =  management 
SUBDIRS += generic
SUBDIRS += importance

SUBLIBS  =  G4biasing_mgt
SUBLIBS +=  G4biasing_gen
SUBLIBS +=  G4biasing_imp

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
