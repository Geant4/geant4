# $Id: GNUmakefile,v 1.18 2010-10-27 07:50:15 gcosmo Exp $
# ----------------------------------------------------------
# GNUmakefile for run library.  Gabriele Cosmo, 3/4/1997.
# ----------------------------------------------------------

name := G4run

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS  = libG4event.lib
GLOBLIBS += libG4tracking.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

ifdef G4FPE_DEBUG
  CPPFLAGS += -DG4FPE_DEBUG
endif

CPPFLAGS += -DG4RUN_ALLOC_EXPORT
CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/processes/scoring/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/transportation/include \
            -I$(G4BASE)/processes/decay/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/digits_hits/detector/include \
            -I$(G4BASE)/digits_hits/hits/include \
            -I$(G4BASE)/digits_hits/digits/include \
            -I$(G4BASE)/digits_hits/utils/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/intercoms/include \
	    -I$(G4BASE)/geometry/biasing/include \
            -I$(G4BASE)/graphics_reps/include 

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

