# --------------------------------------------------------------
# GNUmakefile for examples module.
# --------------------------------------------------------------
name := pythia8_decayer
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all

ifdef PYTHIA8
all: lib bin
  include $(G4INSTALL)/config/binmake.gmk

  INCFLAGS  += -I$(PYTHIA8)/include
  LDLIBS    += -L$(PYTHIA8)/lib -lpythia8
else
all:
	@echo 'ERROR - PYTHIA8 not defined in the environment !'
	@echo '        Tested with Pythia8 release 3.0.5.'
endif
