# $Id: GNUmakefile,v 1.9 2002-06-13 18:01:48 jwellisc Exp $
# ----------------------------------------------------------------
# GNUmakefile for hadronic processes library.  G.Folger 10-Dec-97.
# ----------------------------------------------------------------

name := G4hadronic

SUBDIRS =  management 
SUBDIRS += util
SUBDIRS += processes
SUBDIRS += cross_sections
SUBDIRS += stopping
SUBDIRS += models/low_energy 
SUBDIRS += models/high_energy 
SUBDIRS += models/neutron_hp
SUBDIRS += models/isotope_production
SUBDIRS += models/radiative_decay
SUBDIRS += models/generator/de_excitation
SUBDIRS += models/generator/diffractive_string
SUBDIRS += models/generator/high_energy
SUBDIRS += models/generator/kinetic_model
SUBDIRS += models/generator/management
SUBDIRS += models/generator/pre_equilibrium
SUBDIRS += models/generator/quark_gluon_string
SUBDIRS += models/generator/string_common
SUBDIRS += models/generator/string_fragmentation
SUBDIRS += models/generator/util
SUBDIRS += models/chiral_inv_phase_space/body
SUBDIRS += models/chiral_inv_phase_space/interface
SUBDIRS += models/coherent_elastic
SUBDIRS += models/cascade/evaporation
SUBDIRS += models/cascade/utils
SUBDIRS += models/leading_particle

SUBLIBS =  G4hadronic_mgt G4hadronic_util G4hadronic_proc
SUBLIBS += G4hadronic_xsect G4hadronic_stop  
SUBLIBS += G4hadronic_neu G4hadronic_iso G4hadronic_HE G4hadronic_LE
SUBLIBS += G4hadronic_deex G4hadronic_diffstring G4hadronic_HE_gen G4hadronic_kinetic
SUBLIBS += G4hadronic_man_gen G4hadronic_preequ G4hadronic_radioactivedecay
SUBLIBS += G4hadronic_qgstring G4hadronic_string_common
SUBLIBS += G4hadronic_stringfrag G4hadronic_util_gen
SUBLIBS += G4hadronic_body_ci
SUBLIBS += G4hadronic_interface_ci
SUBLIBS += G4hadronic_coherent_elastic
SUBLIBS += G4hadronic_hetcpp_evaporation
SUBLIBS += G4hadronic_hetcpp_utils
SUBLIBS += G4hadronic_leading_particle

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
