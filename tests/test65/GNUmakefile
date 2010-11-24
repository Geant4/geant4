# $Id: GNUmakefile,v 1.1 2010-11-24 16:15:39 tkoi Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test65
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all

all: lib bin


CPPFLAGS += -I$(G4INSTALL)/source/processes/hadronic/models/lend/include
LDLIBS3 += -lG4had_lend -lexpat

include $(G4INSTALL)/config/binmake.gmk
