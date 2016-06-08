# $Id: GNUmakefile,v 1.3 2002/11/16 17:08:01 hpw Exp $
# ----------------------------------------------------------------
# GNUmakefile for hadronic processes library.  G.Folger 10-Dec-97.
# ----------------------------------------------------------------

name := G4hadronic_gen

SUBDIRS =  de_excitation diffractive_string high_energy scattering kinetic_model management 
SUBDIRS += string_common string_fragmentation util pre_equilibrium quark_gluon_string
SUBDIRS += scattering
SUBLIBS =  G4hadronic_deex G4hadronic_diffstring G4hadronic_HE_gen G4hadronic_kinetic
SUBLIBS += G4hadronic_man_gen G4hadronic_preequ G4hadronic_qgstring G4hadronic_string_common
SUBLIBS += G4hadronic_stringfrag G4hadronic_util_gen G4hadronic_scattering

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/globlib.gmk
