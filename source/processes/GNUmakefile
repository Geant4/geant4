# $Id: GNUmakefile,v 1.11 2001-03-15 09:10:40 gcosmo Exp $
# ----------------------------------------------------------------------
# GNUmakefile for processes library, Gabriele Cosmo 27-Jun-1998
# ----------------------------------------------------------------------

name := G4processes

SUBDIRS  = management optical decay transportation parameterisation
SUBDIRS += electromagnetic/muons electromagnetic/standard
SUBDIRS += electromagnetic/utils electromagnetic/xrays
SUBDIRS += electromagnetic/integral electromagnetic/lowenergy
SUBDIRS += photolepton_hadron
SUBDIRS += hadronic/management hadronic/util hadronic/processes 
SUBDIRS += hadronic/cross_sections hadronic/stopping
SUBDIRS += hadronic/models/low_energy hadronic/models/high_energy
SUBDIRS += hadronic/models/neutron_hp
SUBDIRS += hadronic/models/isotope_production
SUBDIRS += hadronic/models/generator/de_excitation
SUBDIRS += hadronic/models/generator/diffractive_string
SUBDIRS += hadronic/models/generator/high_energy
SUBDIRS += hadronic/models/generator/kinetic_model
SUBDIRS += hadronic/models/generator/management
SUBDIRS += hadronic/models/generator/pre_equilibrium
SUBDIRS += hadronic/models/generator/quark_gluon_string
SUBDIRS += hadronic/models/generator/string_common
SUBDIRS += hadronic/models/generator/string_fragmentation
SUBDIRS += hadronic/models/generator/util
SUBDIRS += hadronic/models/radiative_decay


SUBLIBS  = G4procman G4optical G4decay G4transportation G4parameterisation
SUBLIBS += G4muons G4emstandard G4emutils G4integral G4xrays
SUBLIBS += G4emlowenergy
SUBLIBS += G4photolepton_hadron
SUBLIBS += G4hadronic_mgt G4hadronic_util G4hadronic_proc
SUBLIBS += G4hadronic_xsect G4hadronic_stop G4hadronic_LE G4hadronic_HE
SUBLIBS += G4hadronic_neu G4hadronic_iso
SUBLIBS += G4hadronic_deex G4hadronic_diffstring G4hadronic_HE_gen G4hadronic_kinetic
SUBLIBS += G4hadronic_qgstring
SUBLIBS += G4hadronic_man_gen G4hadronic_preequ G4hadronic_string_common
SUBLIBS += G4hadronic_stringfrag G4hadronic_util_gen
SUBLIBS += G4hadronic_radioactivedecay

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
