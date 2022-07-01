# ------------------------------------------------------------
# GNUmakefile for events library.  Makoto Asai, 5/9/95.
# ------------------------------------------------------------

name := G4event

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS  = libG4tracking.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4EVENT_ALLOC_EXPORT 
##### CPPFLAGS += -DG4_USESMARTSTACK
CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/parameterisation/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
	    -I$(G4BASE)/geometry/biasing/include \
	    -I$(G4BASE)/digits_hits/detector/include \
	    -I$(G4BASE)/digits_hits/hits/include \
	    -I$(G4BASE)/digits_hits/digits/include \
	    -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/graphics_reps/include

ifdef G4_STORE_TRAJECTORY
  CPPFLAGS  += -DG4_STORE_TRAJECTORY
endif

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib
