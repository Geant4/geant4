# ----------------------------------------------------------------

name := G4hadronic_gen

SUBDIRS =  diffraction
SUBDIRS += hadronization
SUBDIRS += management
SUBDIRS += qgsm

SUBLIBS =  G4hadronic_diffstring 
SUBLIBS += G4hadronic_qgstring 
SUBLIBS += G4hadronic_string_common
SUBLIBS += G4hadronic_stringfrag 

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/globlib.gmk
