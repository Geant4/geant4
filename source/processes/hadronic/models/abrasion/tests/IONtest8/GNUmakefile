# $Id: GNUmakefile,v 1.1 2004-05-27 11:38:54 flei Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := IONtest8
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk


ifeq ($(G4SYSTEM),Linux-egcs)
COMPILER := $(shell g++ -v 2>&1 | grep version | grep egcs)
ifeq ($(findstring 1.1,$(COMPILER)),1.1)
$(G4TMP)/$(G4SYSTEM)/$(name)/MLHadronPhysics.o: src/MLHadronPhysics.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/MLHadronPhysics.cc
endif
endif


