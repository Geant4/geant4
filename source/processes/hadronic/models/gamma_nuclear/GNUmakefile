# $Id: GNUmakefile 98814 2016-08-11 09:05:12Z gcosmo $
# -----------------------------------------------------------
# GNUmakefile for hadronic library.  Gabriele Cosmo, 18/9/96.
# -----------------------------------------------------------

name := G4had_gamm_nuclear

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/management/include/ \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/models/cascade/cascade/include \
	    -I$(G4BASE)/processes/hadronic/models/util/include \
	    -I$(G4BASE)/processes/hadronic/models/lend/include \
	    -I$(G4BASE)/processes/hadronic/models/management/include \
	    -I$(G4BASE)/processes/hadronic/models/parton_string/hadronization/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/diffraction/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/management/include \
            -I$(G4BASE)/processes/hadronic/models/theo_high_energy/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/body/include \
            -I$(G4BASE)/processes/hadronic/models/binary_cascade/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/handler/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/management/include \
            -I$(G4BASE)/processes/hadronic/models/pre_equilibrium/exciton_model/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk

