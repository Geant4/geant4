# $Id: GNUmakefile,v 1.9 2005/03/31 16:04:00 duns Exp $
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

# A01PhysicsList use hadronic-lists, so it needs below include command.
#include hadlist.gmk

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

#
# A01HadronPhysics fails to compile in optimized mode on Linux, Darwin, using g++
# so it needs this special rule. Hans-Peter Wellish.
#
ifeq ($(G4SYSTEM),Linux-g++) 
   A01HADRON_NO_OPTIMIZE = 1
endif
ifeq ($(G4SYSTEM),Linux-g++2) 
   A01HADRON_NO_OPTIMIZE = 1
endif
ifeq ($(G4SYSTEM),Linux-g++3) 
   A01HADRON_NO_OPTIMIZE = 1
endif
ifeq ($(G4SYSTEM),Darwin-g++) 
   A01HADRON_NO_OPTIMIZE = 1
endif
ifeq ($(G4SYSTEM),Darwin-g++2) 
   A01HADRON_NO_OPTIMIZE = 1
endif
ifeq ($(G4SYSTEM),Darwin-g++3) 
   A01HADRON_NO_OPTIMIZE = 1
endif

ifdef A01HADRON_NO_OPTIMIZE
$(G4TMP)/$(G4SYSTEM)/$(name)/A01HadronPhysics.o: src/A01HadronPhysics.cc
	@echo "Compiling A01HadronPhysics.cc without optimization..."
	@$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/A01HadronPhysics.cc
endif

ifdef G4ANALYSIS_USE
   CPPFLAGS += `aida-config --include`
   LDLIBS += `aida-config --lib`
endif
