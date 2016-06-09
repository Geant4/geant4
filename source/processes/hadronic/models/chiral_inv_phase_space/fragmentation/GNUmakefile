# --------------------------------------------------------------------
# GNUmakefile for chips cross-sections granular library.  MK 16-Nov-09
# --------------------------------------------------------------------

name := G4hadronic_fragm_ci

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/body/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/fragmentation/include
           
include $(G4INSTALL)/config/common.gmk
