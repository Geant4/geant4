# $Id: GNUmakefile,v 1.13 2004-03-15 11:14:44 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm3
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

#### G4ANALYSIS_USE := true

ifdef G4ANALYSIS_USE
USE_AIDA := true
#### USE_ROOT := true
endif

ifdef USE_AIDA
  CPPFLAGS += -DUSE_AIDA
endif

ifdef USE_ROOT
  CPPFLAGS += -DUSE_ROOT
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef USE_AIDA
  # for the aida-config command see the README file
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
endif

ifdef USE_ROOT
  CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
  LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs) 
endif

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
