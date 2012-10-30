# $Id: GNUmakefile,v 1.27 2010-11-23 13:35:03 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for hadronic models library. G.Folger 10-Dec-97
# --------------------------------------------------------------

name := G4hadronic_models

SUBDIRS  = abrasion
SUBDIRS += binary_cascade
SUBDIRS += cascade/evaporation
SUBDIRS += cascade/utils
SUBDIRS += cascade/cascade
SUBDIRS += chiral_inv_phase_space/body
SUBDIRS += chiral_inv_phase_space/cross_sections
SUBDIRS += chiral_inv_phase_space/fragmentation
SUBDIRS += chiral_inv_phase_space/interface
SUBDIRS += chiral_inv_phase_space/processes
SUBDIRS += coherent_elastic
SUBDIRS += de_excitation/evaporation
SUBDIRS += de_excitation/fermi_breakup
SUBDIRS += de_excitation/fission
SUBDIRS += de_excitation/gem_evaporation
SUBDIRS += de_excitation/handler
SUBDIRS += de_excitation/management
SUBDIRS += de_excitation/multifragmentation
SUBDIRS += de_excitation/photon_evaporation
SUBDIRS += de_excitation/util
SUBDIRS += de_excitation/ablation
SUBDIRS += em_dissociation
SUBDIRS += high_energy
SUBDIRS += im_r_matrix
SUBDIRS += inclxx/utils
SUBDIRS += inclxx/incl_physics
SUBDIRS += inclxx/interface
SUBDIRS += isotope_production
SUBDIRS += lend
SUBDIRS += lll_fission
SUBDIRS += low_energy
SUBDIRS += management
SUBDIRS += neutron_hp
SUBDIRS += lepto_nuclear
SUBDIRS += parton_string/diffraction
SUBDIRS += parton_string/hadronization
SUBDIRS += parton_string/management
SUBDIRS += parton_string/qgsm
SUBDIRS += pre_equilibrium/exciton_model
SUBDIRS += radioactive_decay
SUBDIRS += rpg
SUBDIRS += quasi_elastic
SUBDIRS += qmd
SUBDIRS += theo_high_energy
SUBDIRS += util

SUBLIBS  =  G4hadronic_abrasion
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
SUBLIBS +=  G4hadronic_qgstring
SUBLIBS +=  G4hadronic_quasi_elastic

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/globlib.gmk
