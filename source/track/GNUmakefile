# $Id: GNUmakefile,v 1.1 1999-01-07 16:14:19 gunter Exp $
# ----------------------------------------------------------
# GNUmakefile for track library.  Gabriele Cosmo, 3/4/1997.
# ----------------------------------------------------------

name := G4track

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/processes/management/include \
	    -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

