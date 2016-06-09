# $Id: GNUmakefile,v 1.6.2.1 2010/09/10 14:23:49 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# --------------------------------------------------------------

name := phantom
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
