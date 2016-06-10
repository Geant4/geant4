# $Id: GNUmakefile 84481 2014-10-16 08:53:55Z gcosmo $
# -----------------------------------------------------------
# GNUmakefile for hadronic library.  Gabriele Cosmo, 18/9/96.
# -----------------------------------------------------------

name := G4had_theo_max

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/management/include/ \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
			-I$(G4BASE)/processes/hadronic/models/de_excitation/include \
			-I$(G4BASE)/processes/hadronic/models/management/include \
			-I$(G4BASE)/processes/hadronic/models/util/include \
			-I$(G4BASE)/processes/hadronic/models/im_r_matrix/include \
			-I$(G4BASE)/processes/hadronic/models/pre_equilibrium/include \
			-I$(G4BASE)/processes/hadronic/models/quasi_elastic/include \
			-I$(G4BASE)/processes/hadronic/models/coherent_elastic/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/materials/include
           
include $(G4INSTALL)/config/common.gmk

