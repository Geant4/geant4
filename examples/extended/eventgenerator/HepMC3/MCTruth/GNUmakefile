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
	@echo '        Tested with HepMC release 2.06.03.'
endif
