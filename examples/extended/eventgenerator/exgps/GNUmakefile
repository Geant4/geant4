# $Id: GNUmakefile,v 1.3 2005/07/19 08:35:49 flei Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := exGPS
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

#G4VIS_BUILD_OPENGLX_DRIVER = 0
#G4VIS_USE_OPENGLX = 0
#LDFLAGS += -static

ifdef G4ANALYSIS_USE
   CPPFLAGS += `aida-config --include`
   LDLIBS += `aida-config --lib`
endif
