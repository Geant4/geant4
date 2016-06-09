# $Id: GNUmakefile,v 1.5 2005-06-28 08:37:43 gcosmo Exp $
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

