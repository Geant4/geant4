# --------------------------------------------------------------
# GNUmakefile for chips libraries. MK 19-Nov-09
# --------------------------------------------------------------

name := G4hadronic_chips

SUBDIRS   =  body
SUBDIRS  +=  cross_sections
SUBDIRS  +=  fragmentation
SUBDIRS  +=  processes
SUBDIRS  +=  interface

SUBLIBS   = G4hadronic_body_ci
SUBLIBS  += G4hadronic_crosec_ci
SUBLIBS  += G4hadronic_fragm_ci
SUBLIBS  += G4hadronic_proc_ci
SUBLIBS  += G4hadronic_interface_ci

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif


include $(G4INSTALL)/config/globlib.gmk
