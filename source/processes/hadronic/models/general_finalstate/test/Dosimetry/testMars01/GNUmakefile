# $Id: GNUmakefile,v 1.1 2003-10-08 12:32:20 hpw Exp $
# --------------------------------------------------------------
# GNUmakefile for test module.  Naoko Kanaya, 13/Dec/2001.
# --------------------------------------------------------------

name := testMars01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk
CPPFLAGS  += -I../../include
EXTRALIBS += -lG4hadronic_leading_particle

include $(G4INSTALL)/config/binmake.gmk
