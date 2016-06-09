# $Id: GNUmakefile,v 1.6 2004/06/03 13:55:20 gcosmo Exp $
# -----------------------------------------------------------------------
# GNUmakefile for interfaces.  John Allison, 5/7/98.
# -----------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4interfaces

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4intercoms.lib libG4global.lib

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
