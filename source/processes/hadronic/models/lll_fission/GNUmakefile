# $Id: GNUmakefile 69860 2013-05-16 14:39:02Z gcosmo $
# ------------------------------------------------------------------------
# GNUmakefile for Livermore fission library.  Dennis Wright, 14 March 2007
# ------------------------------------------------------------------------

name := G4had_lll_fis

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/models/neutron_hp/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/processes/hadronic/models/management/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include
           
include $(G4INSTALL)/config/common.gmk

