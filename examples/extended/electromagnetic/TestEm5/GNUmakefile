# $Id: GNUmakefile,v 1.15 2008-06-11 22:15:20 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm5
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: hbook lib bin

#### G4_USE_HBOOK := true
include GNUmakefile.tools_hbook

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/HistoManager.o
			
visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
