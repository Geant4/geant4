# $Id: GNUmakefile,v 1.2 2010-11-03 17:30:47 yarba Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test48
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../
endif

.PHONY: all 

all: lib bin

include $(G4INSTALL)/config/architecture.gmk
#include misc/architecture.gmk

#CPPFLAGS += -I$(G4WORKDIR)/source/particles/management/include
CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
#LDFLAGS  += -L$(G4LIBDIR)
LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)

include $(G4INSTALL)/config/binmake.gmk
#include misc/binmake.gmk

