# $Id: GNUmakefile,v 1.2 2002-07-26 14:53:32 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test31
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

OGLLIBS    := -L$(OGLHOME)/lib -lGLU -lGL
LOADLIBS   += $(OGLLIBS)
CPPFLAGS   += `aida-config --useHBook --include`
LDFLAGS    += `aida-config --useHBook --lib`
#EXTRALIBS += -L$(OGLHOME)/lib -lGL

include $(G4INSTALL)/config/binmake.gmk