# -------------------------------------------------------------
# GNUmakefile for interfaces/common library.  John Allison, 6/7/98.

name := G4UIcommon

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

# For debug mode
# CPPFLAGS += -DG4DEBUG_INTERFACES_COMMON

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

# NO need OGL, then reset OGLGLAGS
OGLFLAGS =
OGLLIBS  =

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/interfaces/basic/include

# Locally adjust source list before including common.gmk
sources := $(addprefix src/,G4InteractorMessenger.cc G4VBasicShell.cc G4VInteractiveSession.cc G4VInteractorManager.cc)

ifdef G4INTY_BUILD_XT
  sources += src/G4Xt.cc
endif

ifdef G4INTY_BUILD_QT
  sources += src/G4Qt.cc
endif

ifdef G4INTY_BUILD_WIN32
  sources += src/G4Win32.cc
endif

include $(G4INSTALL)/config/common.gmk
