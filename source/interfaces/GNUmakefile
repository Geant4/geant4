# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
# -----------------------------------------------------------------------
# GNUmakefile for interfaces.  John Allison, 5/7/98.
# -----------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4interfaces

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4intercoms.lib libG4global.lib

# Definition of macro for moc files (for Qt use only)
MOC_MACRO = -DG4UI_BUILD_QT_SESSION

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

# Libraries which can be used without external packages or libraries...
ifdef G4UI_BUILD
  SUBDIRS := common
  SUBLIBS := G4UIcommon
  SUBDIRS += basic
  SUBLIBS += G4UIbasic
  SUBDIRS += GAG
  SUBLIBS += G4UIGAG
endif

include $(G4INSTALL)/config/globlib.gmk
