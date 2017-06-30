# $Id: GNUmakefile 103669 2017-04-20 15:46:39Z gcosmo $
# ----------------------------------------------------------------------
# GNUmakefile for processes library, Gabriele Cosmo 27-Jun-1998
# ----------------------------------------------------------------------

name := G4processes

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS  = libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib libG4analysis.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib
GLOBLIBS += libG4expat.lib libG4zlib.lib

include $(G4INSTALL)/config/architecture.gmk

SUBDIRS  = management 
SUBDIRS += cuts 
SUBDIRS += solidstate/channeling
SUBDIRS += solidstate/phonon
SUBDIRS += optical 
SUBDIRS += decay 
SUBDIRS += scoring
SUBDIRS += biasing/management
SUBDIRS += biasing/generic
SUBDIRS += biasing/importance
SUBDIRS += transportation 
SUBDIRS += parameterisation
SUBDIRS += electromagnetic/adjoint 
SUBDIRS += electromagnetic/muons 
SUBDIRS += electromagnetic/standard
SUBDIRS += electromagnetic/utils 
SUBDIRS += electromagnetic/xrays
SUBDIRS += electromagnetic/lowenergy
SUBDIRS += electromagnetic/dna/processes
SUBDIRS += electromagnetic/dna/models
SUBDIRS += electromagnetic/dna/utils
SUBDIRS += electromagnetic/dna/management
SUBDIRS += electromagnetic/dna/molecules/management
SUBDIRS += electromagnetic/dna/molecules/types
SUBDIRS += electromagnetic/pii
SUBDIRS += electromagnetic/highenergy
SUBDIRS += electromagnetic/polarisation

SUBLIBS  = G4procman 
SUBLIBS += G4cuts 
SUBLIBS += G4optical 
SUBLIBS += G4channeling
SUBLIBS += G4phonon
SUBLIBS += G4decay
SUBLIBS += G4scoring
SUBLIBS += G4biasing_mgt
SUBLIBS += G4biasing_gen
SUBLIBS += G4biasing_imp
SUBLIBS += G4transportation 
SUBLIBS += G4parameterisation
SUBLIBS += G4muons 
SUBLIBS += G4emstandard 
SUBLIBS += G4emutils 
SUBLIBS += G4xrays 
SUBLIBS += G4emlowenergy
SUBLIBS += G4emhighenergy
SUBLIBS += G4emadjoint
SUBLIBS += G4emdna-processes
SUBLIBS += G4emdna-models
SUBLIBS += G4emdna-utils
SUBLIBS += G4emdna-man
SUBLIBS += G4emdna-molman
SUBLIBS += G4emdna-moltypes
SUBLIBS += G4empii
SUBLIBS += G4empolar

SUBDIRS += hadronic/cross_sections
SUBDIRS += hadronic/management 
SUBDIRS += hadronic/models/abla
SUBDIRS += hadronic/models/abrasion
SUBDIRS += hadronic/models/binary_cascade
SUBDIRS += hadronic/models/cascade/cascade
SUBDIRS += hadronic/models/coherent_elastic
SUBDIRS += hadronic/models/de_excitation/ablation
SUBDIRS += hadronic/models/de_excitation/evaporation
SUBDIRS += hadronic/models/de_excitation/fermi_breakup
SUBDIRS += hadronic/models/de_excitation/fission
SUBDIRS += hadronic/models/de_excitation/gem_evaporation
SUBDIRS += hadronic/models/de_excitation/handler
SUBDIRS += hadronic/models/de_excitation/management
SUBDIRS += hadronic/models/de_excitation/multifragmentation
SUBDIRS += hadronic/models/de_excitation/photon_evaporation
SUBDIRS += hadronic/models/de_excitation/util
SUBDIRS += hadronic/models/em_dissociation
SUBDIRS += hadronic/models/fission
SUBDIRS += hadronic/models/im_r_matrix
SUBDIRS += hadronic/models/inclxx/utils
SUBDIRS += hadronic/models/inclxx/incl_physics
SUBDIRS += hadronic/models/inclxx/interface
SUBDIRS += hadronic/models/lend
SUBDIRS += hadronic/models/management
SUBDIRS += hadronic/models/particle_hp
SUBDIRS += hadronic/models/parton_string/diffraction
SUBDIRS += hadronic/models/parton_string/hadronization
SUBDIRS += hadronic/models/parton_string/management
SUBDIRS += hadronic/models/parton_string/qgsm
SUBDIRS += hadronic/models/lepto_nuclear
SUBDIRS += hadronic/models/pre_equilibrium/exciton_model
SUBDIRS += hadronic/models/qmd
SUBDIRS += hadronic/models/radioactive_decay
SUBDIRS += hadronic/models/rpg
SUBDIRS += hadronic/models/quasi_elastic
SUBDIRS += hadronic/models/theo_high_energy
SUBDIRS += hadronic/models/util
SUBDIRS += hadronic/processes
SUBDIRS += hadronic/stopping
SUBDIRS += hadronic/util

SUBLIBS +=  G4hadronic_abla_abla
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
SUBLIBS +=  G4had_lend
SUBLIBS +=  G4had_fission
SUBLIBS +=  G4had_mod_man
SUBLIBS +=  G4had_part_hp
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

include $(G4INSTALL)/config/globlib.gmk
