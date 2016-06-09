# $Id: GNUmakefile,v 1.31 2010-11-23 13:32:16 gcosmo Exp $
# ----------------------------------------------------------------
# GNUmakefile for hadronic processes library.  G.Folger 10-Dec-97.
# ----------------------------------------------------------------

name := G4hadronic

SUBDIRS =  management 
SUBDIRS += util
SUBDIRS += processes
SUBDIRS += cross_sections
SUBDIRS += stopping


SUBDIRS += models/abrasion
SUBDIRS += models/binary_cascade
SUBDIRS += models/cascade/evaporation
SUBDIRS += models/cascade/utils
SUBDIRS += models/cascade/cascade
SUBDIRS += models/chiral_inv_phase_space/body
SUBDIRS += models/chiral_inv_phase_space/cross_sections
SUBDIRS += models/chiral_inv_phase_space/fragmentation
SUBDIRS += models/chiral_inv_phase_space/interface
SUBDIRS += models/chiral_inv_phase_space/processes
SUBDIRS += models/coherent_elastic
SUBDIRS += models/de_excitation/evaporation
SUBDIRS += models/de_excitation/fermi_breakup
SUBDIRS += models/de_excitation/fission
SUBDIRS += models/de_excitation/gem_evaporation
SUBDIRS += models/de_excitation/handler
SUBDIRS += models/de_excitation/management
SUBDIRS += models/de_excitation/multifragmentation
SUBDIRS += models/de_excitation/photon_evaporation
SUBDIRS += models/de_excitation/util
SUBDIRS += models/de_excitation/ablation
SUBDIRS += models/em_dissociation
SUBDIRS += models/high_energy
SUBDIRS += models/im_r_matrix
SUBDIRS += models/inclxx/utils
SUBDIRS += models/inclxx/incl_physics
SUBDIRS += models/inclxx/interface
SUBDIRS += models/isotope_production
SUBDIRS += models/lend
SUBDIRS += models/lepto_nuclear
SUBDIRS += models/lll_fission
SUBDIRS += models/low_energy
SUBDIRS += models/management
SUBDIRS += models/neutron_hp
SUBDIRS += models/parton_string/diffraction
SUBDIRS += models/parton_string/hadronization
SUBDIRS += models/parton_string/management
SUBDIRS += models/parton_string/qgsm
SUBDIRS += models/pre_equilibrium/exciton_model
SUBDIRS += models/radioactive_decay
SUBDIRS += models/rpg
SUBDIRS += models/quasi_elastic
SUBDIRS += models/qmd
SUBDIRS += models/theo_high_energy
SUBDIRS += models/util

SUBLIBS  =  G4hadronic_abrasion
SUBLIBS +=  G4hadronic_xsect 
SUBLIBS +=  G4hadronic_mgt
SUBLIBS +=  G4hadronic_proc
SUBLIBS +=  G4hadronic_stop
SUBLIBS +=  G4hadronic_util
SUBLIBS +=  G4hadronic_binary
SUBLIBS +=  G4hadronic_coherent_elastic
SUBLIBS +=  G4hadronic_HE
SUBLIBS +=  G4had_im_r_matrix
SUBLIBS +=  G4hadronic_iso
SUBLIBS +=  G4hadronic_LE
SUBLIBS +=  G4hadronic_RPG
SUBLIBS +=  G4hadronic_qmd
SUBLIBS +=  G4had_mod_man
SUBLIBS +=  G4had_neu_hp
SUBLIBS +=  G4had_lend
SUBLIBS +=  G4had_lll_fis
SUBLIBS +=  G4hadronic_radioactivedecay
SUBLIBS +=  G4had_theo_max
SUBLIBS +=  G4had_mod_util
SUBLIBS +=  G4hadronic_bert_cascade
SUBLIBS +=  G4hadronic_hetcpp_evaporation
SUBLIBS +=  G4hadronic_hetcpp_utils
SUBLIBS +=  G4hadronic_inclxx_utils
SUBLIBS +=  G4hadronic_inclxx_physics
SUBLIBS +=  G4hadronic_inclxx_interface
SUBLIBS +=  G4hadronic_body_ci
SUBLIBS +=  G4hadronic_crosec_ci
SUBLIBS +=  G4hadronic_fragm_ci
SUBLIBS +=  G4hadronic_interface_ci
SUBLIBS +=  G4hadronic_proc_ci
SUBLIBS +=  G4hadronic_deex_evaporation
SUBLIBS +=  G4hadronic_deex_fermi_breakup
SUBLIBS +=  G4hadronic_deex_fission
SUBLIBS +=  G4hadronic_deex_gem_evaporation
SUBLIBS +=  G4hadronic_deex_handler
SUBLIBS +=  G4hadronic_deex_management
SUBLIBS +=  G4hadronic_deex_multifragmentation
SUBLIBS +=  G4hadronic_deex_photon_evaporation
SUBLIBS +=  G4hadronic_deex_util
SUBLIBS +=  G4hadronic_ablation
SUBLIBS +=  G4hadronic_em_dissociation
SUBLIBS +=  G4had_preequ_exciton
SUBLIBS +=  G4had_lept_nuclear
SUBLIBS +=  G4had_string_diff
SUBLIBS +=  G4had_string_frag
SUBLIBS +=  G4had_string_man
SUBLIBS +=  G4hadronic_quasi_elastic
SUBLIBS +=  G4hadronic_qgstring

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
