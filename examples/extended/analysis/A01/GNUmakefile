# $Id: GNUmakefile,v 1.5 2003-10-14 20:47:39 duns Exp $
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
ifeq ($(subst g++3,g++,$(subst g++2,g++,$(G4SYSTEM))),Linux-g++)
$(G4TMP)/$(G4SYSTEM)/$(name)/A01HadronPhysics.o: src/A01HadronPhysics.cc
	@echo "Compiling A01HadronPhysics.cc without optimization..."
	@$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/A01HadronPhysics.cc
endif

ifdef G4ANALYSIS_USE
   CPPFLAGS += `aida-config --include`
   LOADLIBS += `aida-config --lib`
endif
