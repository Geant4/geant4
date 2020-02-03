# --------------------------------------------------------------
# GNUmakefile for examples module.
# --------------------------------------------------------------

name := load_gdml
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all
all: lib bin
ifndef XERCESCROOT	
	@echo XERCESCROOT not defined!
endif

include $(G4INSTALL)/config/binmake.gmk

visclean: clean
	@rm -f g4*.prim g4*.eps g4*.wrl
	@rm -f .DAWN_*

clean_all: visclean
	@rm -f g01.gdml
