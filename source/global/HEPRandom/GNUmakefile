# $Id: GNUmakefile 68844 2013-04-08 07:21:18Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for global HEPRandom.  Gunter Folger 12-Dec-2000.
# --------------------------------------------------------------

name := G4heprandom

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GLOB_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/HepGeometry/include
CPPFLAGS += -I$(G4BASE)/global/management/include

include $(G4INSTALL)/config/common.gmk
