# $Id: GNUmakefile,v 1.8 2003-04-09 09:41:09 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test33
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all

all: lib bin

include $(G4INSTALL)/config/architecture.gmk


#LDFLAGS += -lG4hepgeometry

include $(G4INSTALL)/config/binmake.gmk
