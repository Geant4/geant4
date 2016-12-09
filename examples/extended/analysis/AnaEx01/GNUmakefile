# $Id: GNUmakefile 99607 2016-09-28 13:33:42Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := AnaEx01
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

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += -I./shared/include

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

histclean:
	rm ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/HistoManager.o
