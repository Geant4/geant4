# $Id: GNUmakefile,v 1.3 2005-02-22 10:11:09 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for materials library.  Gabriele Cosmo, 18/9/96.
# --------------------------------------------------------------

name := G4materials

GLOBLIBS = libG4global.lib

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/intercoms/include
	    	    
include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

