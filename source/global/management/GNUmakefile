# $Id: GNUmakefile,v 1.3 2001-04-05 09:02:19 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for global management.  Gabriele Cosmo, 26/9/96.
# --------------------------------------------------------------

name := G4globman

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/common.gmk

# Fix for Solaris CC 5.2 Patch 109508-02 to compile without optimization.
#
CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

ifeq ($(G4SYSTEM),SUN-CC5)
COMPILER := $(shell CC -V 2>&1)
ifeq ($(findstring 109508-02,$(COMPILER)),109508-02)
$(G4TMP)/$(G4SYSTEM)/$(name)/G4UnitsTable.o: src/G4UnitsTable.cc
	@echo Compiling G4UnitsTable.cc ...
	@$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4UnitsTable.cc
endif
endif

G4DUMMY_VARIABLE:=$(shell [ ! -d $(G4INCLUDE)/g4std ] && mkdir -p $(G4INCLUDE)/g4std )

installed_includes:=$(foreach file,$(wildcard include/g4std/*),$(shell test -f $(file) && echo $(file)))

installed_includes:=$(patsubst include/%,$(G4INCLUDE)/%,$(installed_includes))

includes::	$(installed_includes)

$(installed_includes): $(G4INCLUDE)/% : include/%
	@cp -p $< $@
