# $Id: GNUmakefile,v 1.1 1999-01-07 16:09:39 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for materials library.  Gabriele Cosmo, 18/9/96.
# --------------------------------------------------------------

name := G4materials

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

