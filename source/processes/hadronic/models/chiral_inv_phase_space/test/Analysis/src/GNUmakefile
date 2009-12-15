# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables for binmake.gmk.
#
INCFLAGS := \
  -I$(G4BASE)/physics_lists/lists/include \
  -I$(G4BASE)/physics_lists/builders/include \
  -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/body/include \
  -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/cross_sections/include \
  -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/fragmentation/include \
  -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/processes/include \
  -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/interface/include \
  -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/test/Analysis/src \
  -I$(G4BASE)/processes/hadronic/models/generator/de_excitation/include \
  -I$(G4BASE)/processes/hadronic/models/generator/pre_equilibrium/include \
  -I$(G4BASE)/processes/hadronic/models/generator/kinetic_model/include \
  -I$(G4BASE)/processes/hadronic/models/de_excitation/util/include \
  -I$(G4BASE)/processes/hadronic/models/de_excitation/fission/include \
  -I$(G4BASE)/processes/hadronic/models/de_excitation/handler/include \
  -I$(G4BASE)/processes/hadronic/models/de_excitation/management/include \
  -I$(G4BASE)/processes/hadronic/models/de_excitation/evaporation/include \
  -I$(G4BASE)/processes/hadronic/models/de_excitation/fermi_breakup/include \
  -I$(G4BASE)/processes/hadronic/models/de_excitation/multifragmentation/include \
  -I$(G4BASE)/processes/hadronic/models/de_excitation/photon_evaporation/include \
  -I$(G4BASE)/processes/hadronic/models/pre_equilibrium/exciton_model/include \
  -I$(G4BASE)/processes/hadronic/models/parton_string/management/include \
  -I$(G4BASE)/processes/hadronic/models/parton_string/hadronization/include \
  -I$(G4BASE)/processes/hadronic/models/parton_string/diffraction/include \
  -I$(G4BASE)/processes/hadronic/models/parton_string/qgsm/include \
  -I$(G4BASE)/processes/hadronic/models/theo_high_energy/include \
  -I$(G4BASE)/processes/hadronic/models/binary_cascade/include \
  -I$(G4BASE)/processes/hadronic/models/management/include \
  -I$(G4BASE)/processes/hadronic/models/util/include \
  -I$(G4BASE)/processes/hadronic/cross_sections/include \
  -I$(G4BASE)/processes/hadronic/management/include \
  -I$(G4BASE)/processes/hadronic/stopping/include \
  -I$(G4BASE)/processes/hadronic/models/low_energy/include \
  -I$(G4BASE)/processes/hadronic/models/high_energy/include \
  -I$(G4BASE)/processes/hadronic/models/generator/quark_gluon_string/include \
  -I$(G4BASE)/processes/hadronic/models/generator/string_common/include \
  -I$(G4BASE)/processes/hadronic/models/generator/string_fragmentation/include \
  -I$(G4BASE)/processes/hadronic/models/generator/diffractive_string/include \
  -I$(G4BASE)/processes/hadronic/models/generator/util/include \
  -I$(G4BASE)/processes/hadronic/models/generator/management/include \
  -I$(G4BASE)/processes/hadronic/models/generator/high_energy/include \
  -I$(G4BASE)/processes/hadronic/processes/include \
  -I$(G4BASE)/processes/hadronic/stopping/include \
  -I$(G4BASE)/processes/hadronic/util/include \
  -I$(G4BASE)/processes/hadronic/models/coherent_elastic/include \
  -I$(G4BASE)/processes/hadronic/models/neutron_hp/include \
  -I$(G4BASE)/processes/hadronic/models/im_r_matrix/include \
  -I$(G4BASE)/processes/hadronic/models/cascade/cascade/include \
  -I$(G4BASE)/processes/hadronic/models/photolepton_hadron/muon_nuclear/include \
  -I$(G4BASE)/processes/electromagnetic/include \
  -I$(G4BASE)/processes/electromagnetic/xrays/include \
  -I$(G4BASE)/processes/decay/include \
  -I$(G4BASE)/processes/cuts/include \
  -I$(G4BASE)/run/include \
  -I$(G4BASE)/event/include \
  -I$(G4BASE)/digits_hits/hits/include \
  -I$(G4BASE)/digits_hits/digits/include \
  -I$(G4BASE)/events+tracks/include \
  -I$(G4BASE)/geometry/management/include \
  -I$(G4BASE)/geometry/navigation/include \
  -I$(G4BASE)/geometry/volumes/include \
  -I$(G4BASE)/geometry/magneticfield/include \
  -I$(G4BASE)/geometry/solids/CSG/include \
  -I$(G4BASE)/global/management/include \
  -I$(G4BASE)/global/HEPRandom/include \
  -I$(G4BASE)/global/HEPGeometry/include \
  -I$(G4BASE)/graphics_reps/include \
  -I$(G4BASE)/intercoms/include \
  -I$(G4BASE)/materials/include \
  -I$(G4BASE)/particles/management/include \
  -I$(G4BASE)/particles/bosons/include \
  -I$(G4BASE)/particles/leptons/include \
  -I$(G4BASE)/particles/hadrons/mesons/include \
  -I$(G4BASE)/particles/hadrons/barions/include \
  -I$(G4BASE)/particles/hadrons/ions/include \
  -I$(G4BASE)/processes/general/include \
  -I$(G4BASE)/processes/management/include \
  -I$(G4BASE)/track/include \
  -I$(G4BASE)/tracking/include 

include $(G4INSTALL)/config/binmake.gmk

