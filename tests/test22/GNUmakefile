# $Id: GNUmakefile,v 1.2 2002-12-09 11:25:10 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test22
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

ifeq ($(G4SYSTEM),HP-aCC)
COMPILER := $(shell aCC -V 2>&1)
ifeq ($(COMPILER), aCC: HP ANSI C++ B3910B A.01.15)
$(G4TMP)/$(G4SYSTEM)/$(name)/Tst02HadronPhysics.o: src/Tst02HadronPhysics.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/Tst02HadronPhysics.cc
endif
endif

ifeq ($(G4SYSTEM),Linux-g++)
COMPILER := $(shell g++ -v 2>&1 | grep version | grep egcs)
ifeq ($(findstring 1.1,$(COMPILER)),1.1)
$(G4TMP)/$(G4SYSTEM)/$(name)/Tst02HadronPhysics.o: src/Tst02HadronPhysics.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/Tst02HadronPhysics.cc
endif
endif
