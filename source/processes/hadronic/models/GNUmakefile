# $Id: GNUmakefile,v 1.10 2002-11-16 17:07:21 hpw Exp $
# --------------------------------------------------------------
# GNUmakefile for hadronic models library. G.Folger 10-Dec-97
# --------------------------------------------------------------

name := G4hadronic_models

SUBDIRS  = low_energy 
SUBDIRS += high_energy 
SUBDIRS += neutron_hp
SUBDIRS += isotope_production
SUBDIRS += generator/de_excitation
SUBDIRS += generator/diffractive_string
SUBDIRS += generator/high_energy
SUBDIRS += generator/kinetic_model
SUBDIRS += generator/management
SUBDIRS += generator/pre_equilibrium
SUBDIRS += generator/quark_gluon_string
SUBDIRS += generator/scattering
SUBDIRS += generator/string_common
SUBDIRS += generator/string_fragmentation
SUBDIRS += generator/util
SUBDIRS += radiative_decay
SUBDIRS += chiral_inv_phase_space/body
SUBDIRS += chiral_inv_phase_space/interface
SUBDIRS += coherent_elastic
SUBDIRS += cascade/evaporation
SUBDIRS += cascade/utils
SUBDIRS += cascade/cascade
SUBDIRS += leading_particle

SUBLIBS  = G4hadronic_HE 
SUBLIBS += G4hadronic_LE 
SUBLIBS += G4hadronic_neu 
SUBLIBS += G4hadronic_iso
SUBLIBS += G4hadronic_deex 
SUBLIBS += G4hadronic_diffstring 
SUBLIBS += G4hadronic_HE_gen 
SUBLIBS += G4hadronic_kinetic
SUBLIBS += G4hadronic_man_gen 
SUBLIBS += G4hadronic_preequ 
SUBLIBS += G4hadronic_radioactivedecay
SUBLIBS += G4hadronic_qgstring 
SUBLIBS += G4hadronic_scattering
SUBLIBS += G4hadronic_string_common 
SUBLIBS += G4hadronic_stringfrag 
SUBLIBS += G4hadronic_util_gen
SUBLIBS += G4hadronic_body_ci
SUBLIBS += G4hadronic_interface_ci
SUBLIBS += G4hadronic_coherent_elastic
SUBLIBS += G4hadronic_hetcpp_evaporation
SUBLIBS += G4hadronic_hetcpp_utils
SUBLIBS += G4hadronic_bert_cascade
SUBLIBS += G4hadronic_leading_particle

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/globlib.gmk
