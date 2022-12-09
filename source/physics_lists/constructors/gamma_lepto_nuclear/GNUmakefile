# ---------------------------------------------------------------------------
# GNUmakefile for physics_lists/constructors/electromagnetic library.  
#  Gunter Folger 10-Jan-2012.
# ---------------------------------------------------------------------------

name := G4phys_ctor_glnuclear

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk


CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/run/include \
	    -I$(G4BASE)/physics_lists/builders/include \
	    -I$(G4BASE)/physics_lists/constructors/electromagnetic/include \
	    -I$(G4BASE)/physics_lists/constructors/factory/include \
	    -I$(G4BASE)/physics_lists/util/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/transportation/include \
            -I$(G4BASE)/processes/decay/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/processes/electromagnetic/highenergy/include \
            -I$(G4BASE)/processes/electromagnetic/muons/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/xrays/include \
            -I$(G4BASE)/processes/optical/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/models/binary_cascade/include \
            -I$(G4BASE)/processes/hadronic/models/cascade/cascade/include \
	    -I$(G4BASE)/processes/hadronic/models/qmd/include \
            -I$(G4BASE)/processes/hadronic/models/coherent_elastic/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/evaporation/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/fermi_breakup/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/handler/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/management/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/multifragmentation/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/photon_evaporation/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/util/include \
	    -I$(G4BASE)/processes/hadronic/models/gamma_nuclear/include \
            -I$(G4BASE)/processes/hadronic/models/lend/include \
            -I$(G4BASE)/processes/hadronic/models/lepto_nuclear/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/diffraction/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/hadronization/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/management/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/qgsm/include \
            -I$(G4BASE)/processes/hadronic/models/pre_equilibrium/exciton_model/include \
            -I$(G4BASE)/processes/hadronic/models/radioactive_decay/include \
            -I$(G4BASE)/processes/hadronic/models/theo_high_energy/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/util/include

include $(G4INSTALL)/config/common.gmk
