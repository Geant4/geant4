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
CPPFLAGS += -I$(G4BASE)/intercoms/include

# Locally adjust source list before including common.gmk
sources := $(addprefix src/,G4UIArrayString.cc G4UIExecutive.cc G4UIcsh.cc G4UIterminal.cc G4VUIshell.cc)

ifeq (,$(findstring WIN32,$(G4SYSTEM)))
  sources += src/G4UItcsh.cc
endif

ifdef G4UI_BUILD_XM_SESSION
  sources += src/G4UIXm.cc
endif

ifdef G4UI_BUILD_QT_SESSION
  sources += src/G4UIQt.cc
endif

ifdef G4UI_BUILD_WIN32_SESSION
  sources += src/G4UIWin32.cc
endif

include $(G4INSTALL)/config/moc.gmk
include $(G4INSTALL)/config/common.gmk

