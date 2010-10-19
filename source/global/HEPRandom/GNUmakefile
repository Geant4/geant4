# $Id: GNUmakefile,v 1.4 2010-10-19 12:25:10 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for global HEPRandom.  Gunter Folger 12-Dec-2000.
# This is needed only for generating "installed" include files.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/HepGeometry/include

include $(G4INSTALL)/config/common.gmk
