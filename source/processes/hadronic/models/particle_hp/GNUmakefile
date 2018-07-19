# $Id: GNUmakefile 76457 2013-11-11 03:01:36Z tkoi $
# -----------------------------------------------------------
# GNUmakefile for hadronic library.  Gabriele Cosmo, 18/9/96.
# -----------------------------------------------------------

name := G4had_part_hp

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT 
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
            -I$(G4BASE)/processes/hadronic/models/de_excitation/photon_evaporation/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/management/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/util/include \
            -I$(G4BASE)/processes/hadronic/models/low_energy/include \
            -I$(G4BASE)/processes/hadronic/models/util/include \
            -I$(G4BASE)/processes/hadronic/models/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/intercoms/include 
           
ifdef G4LIB_BUILD_ZLIB
  CPPFLAGS += -I$(G4BASE)/externals/zlib/include
endif

include $(G4INSTALL)/config/common.gmk

ifdef PHP_AS_HP
   CPPFLAGS += -DPHP_AS_HP
endif
