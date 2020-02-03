# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := AnaEx02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: setup clean_setup all
all: lib bin

setup:
	@echo "Copying files from shared"
	@./shared/scripts/copy_files.sh shared

clean_setup:
	@echo "Removing files copied from shared"
	@./shared/scripts/clean_files.sh shared

# ROOT support
CPPFLAGS += -I$(shell root-config --incdir) 
EXTRALIBS = $(shell root-config --glibs)

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

histclean:
	rm ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/HistoManager.o


