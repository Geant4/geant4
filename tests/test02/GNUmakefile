# $Id: GNUmakefile,v 1.5 1999-12-08 16:30:48 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all tools
 
tools:
	cd ../tools/lib; $(MAKE)
	$(MAKE) all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

ifeq ($(G4SYSTEM),HP-aCC)
COMPILER := $(shell aCC -V 2>&1)
ifeq ($(COMPILER), aCC: HP ANSI C++ B3910B A.01.15)
$(G4TMP)/$(G4SYSTEM)/$(name)/Tst02PhysicsList.o: src/Tst02PhysicsList.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/Tst02PhysicsList.cc
endif
endif

ifeq ($(G4SYSTEM),Linux-g++)
COMPILER := $(shell g++ -v 2>&1 | grep version | grep egcs)
ifeq ($(findstring 1.1,$(COMPILER)),1.1)
$(G4TMP)/$(G4SYSTEM)/$(name)/Tst02PhysicsList.o: src/Tst02PhysicsList.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/Tst02PhysicsList.cc
endif
endif
