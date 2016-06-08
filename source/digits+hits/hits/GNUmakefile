# $Id: GNUmakefile,v 1.1 1999/01/07 16:06:29 gunter Exp $
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
