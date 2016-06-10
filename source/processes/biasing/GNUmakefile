# $Id: GNUmakefile 75716 2013-11-05 15:18:03Z gcosmo $
# ----------------------------------------------------------------
# GNUmakefile for biasing processes library.
# ----------------------------------------------------------------

name := G4biasing

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

SUBDIRS = management generic importance

SUBLIBS = G4biasing_mgt G4biasing_gen G4biasing_imp

include $(G4INSTALL)/config/globlib.gmk
