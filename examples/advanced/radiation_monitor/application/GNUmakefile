# $Id: GNUmakefile,v 1.7.4.1 2009/08/11 14:20:35 gcosmo Exp $
# --------------------------------------------------------------------
# GNUmakefile for exaple subpackage
# R. Capra, 09 Sep 2005
# --------------------------------------------------------------------

name := Radmon

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_USE.gmk
include $(G4INSTALL)/config/G4VIS_USE.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS +=  -I$(G4BASE)/global/management/include \
              -I$(G4BASE)/global/HEPRandom/include \
              -I$(G4BASE)/global/HEPGeometry/include \
              -I$(G4BASE)/global/HEPNumerics/include \
              -I$(G4BASE)/run/include \
              -I$(G4BASE)/materials/include \
              -I$(G4BASE)/parameterisations/trd_clusters/include \
              -I$(G4BASE)/parameterisations/gflash/include \
              -I$(G4BASE)/particles/management/include \
              -I$(G4BASE)/particles/bosons/include \
              -I$(G4BASE)/particles/leptons/include \
              -I$(G4BASE)/particles/hadrons/barions/include \
              -I$(G4BASE)/particles/hadrons/ions/include \
              -I$(G4BASE)/particles/hadrons/mesons/include \
              -I$(G4BASE)/particles/shortlived/include \
              -I$(G4BASE)/processes/management/include \
              -I$(G4BASE)/processes/cuts/include \
              -I$(G4BASE)/processes/decay/include \
              -I$(G4BASE)/processes/optical/include \
              -I$(G4BASE)/processes/transportation/include \
              -I$(G4BASE)/processes/parameterisation/include \
              -I$(G4BASE)/processes/electromagnetic/standard/include \
              -I$(G4BASE)/processes/electromagnetic/muons/include \
              -I$(G4BASE)/processes/electromagnetic/utils/include \
              -I$(G4BASE)/processes/electromagnetic/xrays/include \
              -I$(G4BASE)/processes/electromagnetic/lowenergy/include \
              -I$(G4BASE)/processes/electromagnetic/highenergy/include \
              -I$(G4BASE)/processes/hadronic/cross_sections/include \
	      -I$(G4BASE)/processes/hadronic/management/include \
	      -I$(G4BASE)/processes/hadronic/models/abrasion/include \
              -I$(G4BASE)/processes/hadronic/models/binary_cascade/include \
	      -I$(G4BASE)/processes/hadronic/models/cascade/cascade/include \
	      -I$(G4BASE)/processes/hadronic/models/cascade/evaporation/include \
	      -I$(G4BASE)/processes/hadronic/models/cascade/utils/include \
              -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/body/include \
              -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/interface/include \
              -I$(G4BASE)/processes/hadronic/models/coherent_elastic/include \
	      -I$(G4BASE)/processes/hadronic/models/de_excitation/ablation/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/evaporation/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/fermi_breakup/include \
	      -I$(G4BASE)/processes/hadronic/models/de_excitation/fission/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/gem_evaporation/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/handler/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/management/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/multifragmentation/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/photon_evaporation/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/util/include \
              -I$(G4BASE)/processes/hadronic/models/em_dissociation/include \
              -I$(G4BASE)/processes/hadronic/models/general_finalstate/include \
              -I$(G4BASE)/processes/hadronic/models/high_energy/include \
              -I$(G4BASE)/processes/hadronic/models/im_r_matrix/include \
              -I$(G4BASE)/processes/hadronic/models/inucl/include \
              -I$(G4BASE)/processes/hadronic/models/isotope_production/include \
              -I$(G4BASE)/processes/hadronic/models/leading_particle/include \
              -I$(G4BASE)/processes/hadronic/models/low_energy/include \
              -I$(G4BASE)/processes/hadronic/models/management/include \
              -I$(G4BASE)/processes/hadronic/models/neutron_hp/include \
              -I$(G4BASE)/processes/hadronic/models/parton_string/diffraction/include \
              -I$(G4BASE)/processes/hadronic/models/parton_string/hadronization/include \
              -I$(G4BASE)/processes/hadronic/models/parton_string/management/include \
              -I$(G4BASE)/processes/hadronic/models/parton_string/qgsm/include \
              -I$(G4BASE)/processes/hadronic/models/photolepton_hadron/muon_nuclear/include \
              -I$(G4BASE)/processes/hadronic/models/pre_equilibrium/exciton_model/include \
              -I$(G4BASE)/processes/hadronic/models/quark_molecular_dynamics/body/include \
              -I$(G4BASE)/processes/hadronic/models/quark_molecular_dynamics/util/include \
              -I$(G4BASE)/processes/hadronic/models/radioactive_decay/include \
              -I$(G4BASE)/processes/hadronic/models/relativistic_qmd/include \
              -I$(G4BASE)/processes/hadronic/models/relativistic_qmd/uplusplus/include \
              -I$(G4BASE)/processes/hadronic/models/theo_high_energy/include \
              -I$(G4BASE)/processes/hadronic/models/util/include \
              -I$(G4BASE)/processes/hadronic/processes/include \
              -I$(G4BASE)/processes/hadronic/stopping/include \
              -I$(G4BASE)/processes/hadronic/util/include \
              -I$(G4BASE)/geometry/management/include \
              -I$(G4BASE)/geometry/solids/CSG/include \
              -I$(G4BASE)/geometry/solids/specific/include \
              -I$(G4BASE)/geometry/solids/BREPS/include \
              -I$(G4BASE)/geometry/solids/Boolean/include \
              -I$(G4BASE)/geometry/divisions/include \
              -I$(G4BASE)/geometry/volumes/include \
              -I$(G4BASE)/geometry/navigation/include \
              -I$(G4BASE)/geometry/magneticfield/include \
              -I$(G4BASE)/geometry/biasing/include \
              -I$(G4BASE)/track/include \
              -I$(G4BASE)/tracking/include \
              -I$(G4BASE)/digits_hits/detector/include \
              -I$(G4BASE)/digits_hits/hits/include \
              -I$(G4BASE)/digits_hits/digits/include \
              -I$(G4BASE)/digits_hits/utils/include \
              -I$(G4BASE)/persistency/include \
              -I$(G4BASE)/readout/include \
              -I$(G4BASE)/event/include \
              -I$(G4BASE)/graphics_reps/include \
              -I$(G4BASE)/intercoms/include \
            -I../generalPurpose/include \
            -I../generator/include \
            -I../detector/include \
            -I../materials/include \
            -I../physics/include \
            -I../analysis/include \
            -I../sensitivity/include \
            -I../userActions/include


include $(G4INSTALL)/config/analysis.gmk
include $(G4INSTALL)/config/common.gmk
