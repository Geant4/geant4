# $Id: GNUmakefile,v 1.3 2003/05/28 08:03:22 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for config scripts.  Sergei Sadilov, 13/12/02.
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
