# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
# ----------------------------------------------------------
# GNUmakefile for tracking library.  Katsuya Amako, 5/9/95.
# ----------------------------------------------------------

name := G4error_propagation

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS  = libG4run.lib libG4event.lib
GLOBLIBS += libG4tracking.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/run/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/standard/include \
            -I$(G4BASE)/processes/electromagnetic/muons/include \
            -I$(G4BASE)/processes/transportation/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/digits_hits/hits/include \
            -I$(G4BASE)/digits_hits/digits/include 

ifdef G4EVERBOSE
  CPPFLAGS += -DG4EVERBOSE
endif

include $(G4INSTALL)/config/common.gmk

.PHONY: global

global: lib
