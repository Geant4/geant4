# $Id: GNUmakefile,v 1.1 2002-11-19 10:32:59 murakami Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
name := hepmcEx01-clhep
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all
all : lib bin

include $(G4INSTALL)/config/binmake.gmk

# check gcc
USE_GCC := 0
ifeq ($(G4SYSTEM),Linux-g++)
  USE_GCC := 1
endif
ifeq ($(G4SYSTEM),Linux-egcs)
  USE_GCC := 1
endif

# PYTHIA
  CERNLIB_DIR := /cern/pro/lib

ifeq ($(USE_GCC),1)
  LDLIBS1     += -L$(CERNLIB_DIR) -lpythia6 -lg2c
else
  LDLIBS1     += -L$(CERNLIB_DIR) -lpythia6 # and something?
endif

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

