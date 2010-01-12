# $Id: GNUmakefile,v 1.3 2010-01-12 09:23:04 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.
# --------------------------------------------------------------

name := test37
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk


include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_* *.out


