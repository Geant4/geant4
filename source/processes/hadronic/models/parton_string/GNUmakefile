# --------------------------------------------------------------
# GNUmakefile for chips library. HPW 20-Nov-99
# --------------------------------------------------------------

name := G4had_string


ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

SUBDIRS  = diffraction
SUBDIRS += hadronization   
SUBDIRS += management
SUBDIRS += qgsm


SUBLIBS  = G4had_string_diff
SUBLIBS += G4had_string_frag
SUBLIBS += G4had_string_man
SUBLIBS += G4hadronic_qgstring

include $(G4INSTALL)/config/globlib.gmk
