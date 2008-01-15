# $Id: GNUmakefile,v 1.5 2008-01-15 11:04:26 lgarnier Exp $
# -------------------------------------------------------------
# GNUmakefile for interfaces/basic library.  John Allison, 6/7/98.

name := G4UIbasic

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

# Definition of macro for moc files (for Qt use only)
MOC_MACRO = -DG4UI_BUILD_QT_SESSION

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/interfaces/common/include
CPPFLAGS += -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk

