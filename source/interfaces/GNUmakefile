# $Id: GNUmakefile,v 1.3 2003-10-22 15:17:53 gbarrand Exp $
# -----------------------------------------------------------------------
# GNUmakefile for interfaces.  John Allison, 5/7/98.
# -----------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4interfaces

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

G4LIBDIR := $(G4LIB)/$(G4SYSTEM)
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(name)

# Libraries which can be used without external packages or libraries...
ifdef G4UI_BUILD
  SUBDIRS := common
  SUBLIBS := G4UIcommon
  SUBDIRS += basic
  SUBLIBS += G4UIbasic
  SUBDIRS += GAG
  SUBLIBS += G4UIGAG

########################### interfaces/XVT - not currently implemented #####

endif #G4UI_BUILD

include $(G4INSTALL)/config/globlib.gmk
