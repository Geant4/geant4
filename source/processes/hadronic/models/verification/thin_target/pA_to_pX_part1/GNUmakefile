# $Id: GNUmakefile,v 1.3 2003-05-27 14:35:04 hpw Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := ppinc
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../../../..
endif
  INCFLAGS += -I$(G4BASE)/source/processes/hadronic/models/verification/thin_target/pA_to_pX_part1/include/AIDA

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

