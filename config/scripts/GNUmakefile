# $Id: GNUmakefile,v 1.2 2003/02/21 15:23:26 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := Marker
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all

all: lib bin

include $(G4INSTALL)/config/scripts/binmake.gmk

CXXFLAGS_FINAL := $(CXXFLAGS) 
INCFLAGS_FINAL := $(INCFLAGS)
LDLIBS_FINAL := $(LDLIBS)
LDLIB3_FINAL := $(LDLIBS3)
LDLIB4_FINAL := $(LDLIBS4)
#LDFLAGS_FINAL := $(LDFLAGS)
