# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleB3aExt
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: setup clean_setup all
all: lib bin

setup:
	@echo "Copying files from scoreWriter"
	@cp -rp scoreWriter/include/* include
	@cp -rp scoreWriter/src/* src

clean_setup:
	@echo "Removing files copied from scoreWriter"
	@rm -fr include/G4ScoreNtupleWriter*
	@rm -fr src/G4ScoreNtupleWriter*

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += -I./scoreWriter/include

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
