# $Id: GNUmakefile,v 1.3 2004-03-10 12:32:48 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm7
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

#### USE_AIDA := true
ifndef USE_AIDA
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
