# $Id: GNUmakefile,v 1.2 2002-07-11 08:12:22 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleB01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin
include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4INSTALL)/examples/extended/biasing/B09/include

LDFLAGS += -L$(G4WORKDIR)/tmp/Linux-g++/exampleB09\
	   -lexampleB09

include $(G4INSTALL)/config/binmake.gmk
