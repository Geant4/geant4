# $Id: GNUmakefile,v 1.1 1999-01-07 16:09:20 gunter Exp $
# -------------------------------------------------------------
# GNUmakefile for intercoms library.  Gabriele Cosmo, 11/11/96.
# -------------------------------------------------------------

name := G4intercoms

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

