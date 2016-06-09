# $Id: GNUmakefile,v 1.4 2006/07/05 11:53:02 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
name := hepmcEx01
G4TARGET := $(name)
G4EXLIB := true

include hadronic_lists.gmk

.PHONY: all
all : had_list lib bin

had_list:
	@echo "This example uses the QGSP 'educated guess physics list'"
	@echo "(see $(G4INSTALL)/physics_lists/hadronic)."
	@echo "The hadronic physics lists should be compiled before"
	@echo "compiling this example."

include $(G4INSTALL)/config/binmake.gmk

# -----------------------------------------------------------------
# HepMC and PYTHIA

# if you do not use Pythia library, comment out the next line.
#
G4LIB_USE_PYTHIA := 1
ifdef G4LIB_USE_PYTHIA
  CPPFLAGS += -DG4LIB_USE_PYTHIA
endif

# Define HepMC installation path (HepMC 1.27)
#
ifndef HEPMC_DIR
	@echo HEPMC_DIR not defined!
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
