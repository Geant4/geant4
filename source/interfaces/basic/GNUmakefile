# $Id: GNUmakefile,v 1.7 2009/10/19 07:37:48 lgarnier Exp $
# -------------------------------------------------------------
# GNUmakefile for interfaces/basic library.  John Allison, 6/7/98.

name := G4UIbasic

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

# Definition of macro for moc files (for Qt use only)
MOC_MACRO = -DG4UI_BUILD_QT_SESSION

# For debug mode
# CPPFLAGS += -DG4DEBUG_INTERFACES_BASIC

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/interfaces/common/include
CPPFLAGS += -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/moc.gmk
include $(G4INSTALL)/config/common.gmk

