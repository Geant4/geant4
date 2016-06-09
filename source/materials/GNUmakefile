# $Id: GNUmakefile,v 1.2 2004/01/30 09:57:20 gbarrand Exp $
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
            -I$(G4BASE)/global/HEPGeometry/include 
	    
include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

