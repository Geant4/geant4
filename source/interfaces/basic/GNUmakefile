# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
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

# NO need OGL, then reset OGLGLAGS
OGLFLAGS =
OGLLIBS  =

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/interfaces/common/include
CPPFLAGS += -I$(G4BASE)/interfaces/GAG/include
CPPFLAGS += -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/moc.gmk
include $(G4INSTALL)/config/common.gmk

