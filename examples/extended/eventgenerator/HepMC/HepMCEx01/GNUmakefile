# $Id: GNUmakefile,v 1.2 2002-05-28 14:05:16 murakami Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
name := hepmcEx01
G4TARGET := $(name)
G4EXLIB := true

#ifndef G4WORKDIR
#  G4WORKDIR := .
#endif

.PHONY: all
all : lib bin

include $(G4INSTALL)/config/binmake.gmk

# -----------------------------------------------------------------
# HepMC and PYTHIA

# if you donnot use Pythia library, comment out.
USE_PYTHIA := 1

# HepMC
HEPMC_DIR := /opt/heplib/HepMC/1.21
INCFLAGS  += -I$(HEPMC_DIR)

# check gcc
USE_GCC := 0
ifeq ($(G4SYSTEM),Linux-g++)
  USE_GCC := 1
endif
ifeq ($(G4SYSTEM),Linux-egcs)
  USE_GCC := 1
endif

ifeq ($(USE_GCC),1)
   LDLIBS1  += -Wl,-rpath $(HEPMC_DIR)
endif

ifdef USE_PYTHIA
    LDLIBS1  += -L$(HEPMC_DIR) -lHepMC
else
    LDLIBS1  += -L$(HEPMC_DIR) -lHepMC $(G4TMPDIR)/HEPEvtcom.o
endif

ifdef USE_PYTHIA 
  CERNLIB_DIR := /cern/pro/lib
  LDLIBS1     += -L$(CERNLIB_DIR) -lpythia6

  # if you get pythia initialization errors, try the next flag.
  #NEED_INITPYDATA := 1
  ifdef NEED_INITPYDATA
    CXXFLAGS  += -DNEED_INITPYDATA
    LDLIBS1   += $(HepMC_DIR)/src/initpydata.o
   #LDLIBS1   += external/initpydata.o
  endif

  LDLIBS1     += $(FCLIBS)
endif
# -----------------------------------------------------------------

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

# On HP-acc optimization has to be turned off to compile ExN04HadronPhysics.
#
ifeq ($(G4SYSTEM),HP-aCC)
  COMPILER := $(shell aCC -V 2>&1)
ifeq ($(COMPILER), aCC: HP ANSI C++ B3910B A.01.15)
$(G4TMP)/$(G4SYSTEM)/$(name)/XXHadronPhysics.o: src/XXHadronPhysics.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/XXHadronPhysics.cc
endif
endif

