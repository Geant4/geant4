# $Id: GNUmakefile 67044 2013-01-30 08:50:06Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for materials library.  Gabriele Cosmo, 18/9/96.
# --------------------------------------------------------------

name := G4materials

GLOBLIBS = libG4intercoms.lib libG4global.lib

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

