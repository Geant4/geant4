# $Id: GNUmakefile,v 1.1 1999-01-07 16:09:29 gunter Exp $
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
ifndef G4UI_BUILD
  G4UI_BUILD := yes
endif
ifndef G4UI_BUILD_BASIC_SESSION
  G4UI_BUILD_BASIC_SESSION := yes
endif

########################### interfaces/common #######################

ifdef G4UI_BUILD
  SUBDIRS := common
  SUBLIBS := G4UIcommon
endif

########################### interfaces/basic #######################

ifdef G4UI_BUILD_BASIC_SESSION
  SUBDIRS += basic
  SUBLIBS += G4UIbasic
endif

########################### interfaces/GAG #######################

ifdef G4UI_BUILD_GAG_SESSION
  SUBDIRS += GAG
  SUBLIBS += G4UIGAG
endif

########################### interfaces/OPACS #######################

ifdef G4UI_BUILD_WO_SESSION
  SUBDIRS += OPACS
  SUBLIBS += G4UIOPACS
endif

########################### interfaces/XVT - not currently implemented #####

include $(G4INSTALL)/config/globlib.gmk
