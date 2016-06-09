# --------------------------------------------------------------
# GNUmakefile for chips library. HPW 20-Nov-99
# --------------------------------------------------------------

name := G4hadronic_deex


ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

SUBDIRS  = handler
SUBDIRS += management       
SUBDIRS += evaporation
SUBDIRS += gem_evaporation
SUBDIRS += photon_evaporation
SUBDIRS += fission
SUBDIRS += fermi_breakup
SUBDIRS += multifragmentation
SUBDIRS += util
SUBDIRS += ablation


SUBLIBS  = G4hadronic_deex_handler
SUBLIBS += G4hadronic_deex_management
SUBLIBS += G4hadronic_deex_evaporation
SUBLIBS += G4hadronic_deex_gem_evaporation
SUBLIBS += G4hadronic_deex_photon_evaporation
SUBLIBS += G4hadronic_deex_fission
SUBLIBS += G4hadronic_deex_fermi_breakup
SUBLIBS += G4hadronic_deex_multifagmentation
SUBLIBS += G4hadronic_deex_util
SUBLIBS += G4hadronic_ablation

include $(G4INSTALL)/config/globlib.gmk
