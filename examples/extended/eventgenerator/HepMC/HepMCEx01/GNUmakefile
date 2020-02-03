# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
name := HepMCEx01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all

ifdef HEPMC_DIR
all : pythia6 lib bin
# if you do not use Pythia library, replace the line above with the one below.
# all : lib bin

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
      LDLIBS1  += -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio -L$(G4TMPDIR) -lPythia6 -lg2c
  else
      LDLIBS1  += -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio $(G4TMPDIR)/HEPEvtcom.o
  endif

  # Path for PYTHIA Fortran library. Based on CERNLIB-2005.
  # Add "/cern/pro/bin" to $PATH first !
  #
  #ifdef G4LIB_USE_PYTHIA 
  #  LDLIBS1     += $(shell cernlib -v pro pythia6205 pdflib804 mathlib) -lg2c
  #endif

FCFLAGS += -c

pythia6: $(G4TMPDIR)/libPythia6.so

$(G4TMPDIR)/libPythia6.so: $(G4TMPDIR)/pythia6.o
	$(FC) -shared -Wl,-soname,libPythia6.so -o $(G4TMPDIR)/libPythia6.so  $(G4TMPDIR)/pythia6.o
 
$(G4TMPDIR)/pythia6.o:
	$(FC) $(FCFLAGS) $(PYTHIA6)/pythia-$(PYTHIA6_VERSION).f -o $(G4TMPDIR)/pythia6.o

else

all:
	@echo 'ERROR - HEPMC_DIR not defined in the environment !'
	@echo '        Requires HepMC release 1.27.'
endif
