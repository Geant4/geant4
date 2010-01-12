# $Id: GNUmakefile,v 1.5 2010-01-12 09:45:26 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.
# --------------------------------------------------------------

name := test35
G4TARGET := test35
G4EXLIB  := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all tools

all: lib bin

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Histo.o


