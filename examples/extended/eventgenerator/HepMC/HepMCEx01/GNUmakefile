# $Id: GNUmakefile,v 1.6 2006-11-22 15:09:21 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
name := hepmcEx01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all

ifdef HEPMC_DIR
all : lib bin
  include $(G4INSTALL)/config/binmake.gmk

  # -----------------------------------------------------------------
  # HepMC and PYTHIA

  # if you do not use Pythia library, comment out the next line.
  #
  G4LIB_USE_PYTHIA := 1
  ifdef G4LIB_USE_PYTHIA
    CPPFLAGS += -DG4LIB_USE_PYTHIA
  endif

  INCFLAGS  += -I$(HEPMC_DIR)/include

  ifdef G4LIB_USE_PYTHIA
      LDLIBS1  += -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio
  else
      LDLIBS1  += -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio $(G4TMPDIR)/HEPEvtcom.o
  endif

  # Path for PYTHIA Fortran library. Based on CERNLIB-2005.
  # Add "/cern/pro/bin" to $PATH first !
  #
  ifdef G4LIB_USE_PYTHIA 
    LDLIBS1     += $(shell cernlib -v pro pythia6205 pdflib804 mathlib) -lg2c
  endif
else
all:
	@echo 'ERROR - HEPMC_DIR not defined in the environment !'
	@echo '        Requires HepMC release 1.27.'
endif
