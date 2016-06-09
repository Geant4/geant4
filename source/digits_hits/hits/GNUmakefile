# $Id: GNUmakefile,v 1.1 2003/10/03 10:16:22 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for digits+hits library.  Makoto Asai, 1/11/96.
# --------------------------------------------------------------

name := G4hits

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPGeometry/include

include $(G4INSTALL)/config/common.gmk
