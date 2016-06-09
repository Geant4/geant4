# $Id: GNUmakefile,v 1.8 2010-10-27 07:53:29 gcosmo Exp $
# ----------------------------------------------------------
# GNUmakefile for track library.  Gabriele Cosmo, 3/4/1997.
# ----------------------------------------------------------

name := G4track

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS  = libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4TRACK_ALLOC_EXPORT
CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/management/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

