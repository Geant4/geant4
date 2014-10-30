# $Id: GNUmakefile 66356 2012-12-18 09:02:32Z gcosmo $
# --------------------------------------------------------------------
# GNUmakefile for externals/usolids library.
# --------------------------------------------------------------------

name := G4geomUSolids

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

GLOBLIBS  = 

include $(G4INSTALL)/config/architecture.gmk

ifdef G4USOLDEBUG
  CPPFLAGS += -DG4USOLDEBUG
endif

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib
