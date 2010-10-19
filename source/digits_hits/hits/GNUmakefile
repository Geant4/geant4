# $Id: GNUmakefile,v 1.3 2010-10-19 09:54:03 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for digits+hits library.  Makoto Asai, 1/11/96.
# --------------------------------------------------------------

name := G4hits

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4ALLOC_EXPORT
CPPFLAGS += \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPGeometry/include

include $(G4INSTALL)/config/common.gmk
