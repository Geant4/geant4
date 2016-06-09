# $Id: GNUmakefile,v 1.1 2006/11/22 14:51:26 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# --------------------------------------------------------------
name := mctruthex
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all

ifdef HEPMC_DIR
all: lib bin
  include $(G4INSTALL)/config/binmake.gmk

  INCFLAGS  += -I$(HEPMC_DIR)/include
  LDLIBS    += -L$(HEPMC_DIR)/lib -lHepMC
else
all:
	@echo 'ERROR - HEPMC_DIR not defined in the environment !'
	@echo '        Requires HepMC release 1.27.'
endif
