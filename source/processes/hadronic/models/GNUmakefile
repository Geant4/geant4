# $Id: GNUmakefile,v 1.7 2001-11-15 17:46:55 hpw Exp $
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
SUBDIRS += generator/string_common
SUBDIRS += generator/string_fragmentation
SUBDIRS += generator/util
SUBDIRS += radiative_decay
SUBDIRS += chiral_inv_phase_space/body
SUBDIRS += chiral_inv_phase_space/interface
SUBDIRS += coherent_elastic

SUBLIBS  = G4hadronic_HE G4hadronic_LE G4hadronic_neu G4hadronic_iso
SUBLIBS += G4hadronic_deex G4hadronic_diffstring G4hadronic_HE_gen G4hadronic_kinetic
SUBLIBS += G4hadronic_man_gen G4hadronic_preequ G4hadronic_radioactivedecay
SUBLIBS += G4hadronic_qgstring G4hadronic_string_common
SUBLIBS += G4hadronic_stringfrag G4hadronic_util_gen
SUBLIBS += G4hadronic_body_ci
SUBLIBS += G4hadronic_interface_ci
SUBLIBS += G4hadronic_coherent_elastic

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/globlib.gmk
