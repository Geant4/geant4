# $Id: GNUmakefile 106468 2017-10-11 08:11:47Z gcosmo $
# ----------------------------------------------------------------
# GNUmakefile for hadronic processes library.  G.Folger 10-Dec-97.
# ----------------------------------------------------------------

name := G4hadronic

SUBDIRS =  management 
SUBDIRS += util
SUBDIRS += processes
SUBDIRS += cross_sections
SUBDIRS += stopping
SUBDIRS += models/abla
SUBDIRS += models/abrasion
SUBDIRS += models/binary_cascade
SUBDIRS += models/cascade/cascade
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
SUBDIRS += models/fission
SUBDIRS += models/gamma_nuclear
SUBDIRS += models/im_r_matrix
SUBDIRS += models/inclxx/utils
SUBDIRS += models/inclxx/incl_physics
SUBDIRS += models/inclxx/interface
SUBDIRS += models/lend
SUBDIRS += models/lepto_nuclear
SUBDIRS += models/management
SUBDIRS += models/particle_hp
SUBDIRS += models/parton_string/diffraction
SUBDIRS += models/parton_string/hadronization
SUBDIRS += models/parton_string/management
SUBDIRS += models/parton_string/qgsm
SUBDIRS += models/pre_equilibrium/exciton_model
SUBDIRS += models/radioactive_decay
SUBDIRS += models/rpg
SUBDIRS += models/qmd
SUBDIRS += models/quasi_elastic
SUBDIRS += models/theo_high_energy
SUBDIRS += models/util

SUBLIBS  =  G4hadronic_abla_abla
SUBLIBS +=  G4hadronic_abrasion
SUBLIBS +=  G4hadronic_xsect 
SUBLIBS +=  G4hadronic_mgt
SUBLIBS +=  G4hadronic_proc
SUBLIBS +=  G4hadronic_stop
SUBLIBS +=  G4hadronic_util
SUBLIBS +=  G4hadronic_binary
SUBLIBS +=  G4hadronic_coherent_elastic
SUBLIBS +=  G4had_im_r_matrix
SUBLIBS +=  G4hadronic_RPG
SUBLIBS +=  G4hadronic_qmd
SUBLIBS +=  G4had_mod_man
SUBLIBS +=  G4had_part_hp
SUBLIBS +=  G4had_lend
SUBLIBS +=  G4had_fission
SUBLIBS +=  G4had_gamm_nuclear
SUBLIBS +=  G4hadronic_radioactivedecay
SUBLIBS +=  G4had_theo_max
SUBLIBS +=  G4had_mod_util
SUBLIBS +=  G4hadronic_bert_cascade
SUBLIBS +=  G4hadronic_inclxx_utils
SUBLIBS +=  G4hadronic_inclxx_physics
SUBLIBS +=  G4hadronic_inclxx_interface
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
