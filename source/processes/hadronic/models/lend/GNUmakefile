# $Id: GNUmakefile 94435 2015-11-17 14:24:19Z gcosmo $
# -----------------------------------------------------------
# GNUmakefile for hadronic lend library.
# -----------------------------------------------------------

name := G4had_lend

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4PROCESSES_EXPORT
CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT

ifeq ($(G4SYSTEM), WIN32-VC)
CPPFLAGS += -I$(G4BASE)/externals/expat/include
endif
ifdef G4LIB_BUILD_EXPAT
CPPFLAGS += -I$(G4BASE)/externals/expat/include
endif

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/management/include/ \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/handler/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/photon_evaporation/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/management/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/util/include \
            -I$(G4BASE)/processes/hadronic/models/pre_equilibrium/exciton_model/include \
            -I$(G4BASE)/processes/hadronic/models/low_energy/include \
            -I$(G4BASE)/processes/hadronic/models/util/include \
            -I$(G4BASE)/processes/hadronic/models/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
