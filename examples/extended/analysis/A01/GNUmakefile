# $Id: GNUmakefile,v 1.3 2003-10-13 20:52:06 duns Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := A01app
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

#
# A01HadronPhysics fails to compile in optimized mode on Linux,
# so it needs this special rule. Hans-Peter Wellish.
#
ifeq ($(G4SYSTEM),Linux-g++)
$(G4TMP)/$(G4SYSTEM)/$(name)/A01HadronPhysics.o: src/A01HadronPhysics.cc
	@echo "Compiling A01HadronPhysics.cc ..."
	@$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/A01HadronPhysics.cc
endif
