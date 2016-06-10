# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
# --------------------------------------------------------------
# GNUmakefile for  processes/decay library.  G.Folger 9/12/97
# --------------------------------------------------------------

name := G4cuts

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/materials/include \
	    -I$(G4BASE)/geometry/management/include \
	    -I$(G4BASE)/processes/management/include \
	    -I$(G4BASE)/particles/management/include \
	    -I$(G4BASE)/intercoms/include 

include $(G4INSTALL)/config/common.gmk
