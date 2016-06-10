# $Id: GNUmakefile 85989 2014-11-06 15:40:05Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := AnaEx01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

# HBOOK support
#
#### G4_USE_HBOOK := true

.PHONY: setup clean_setup all
all: hbook lib bin

setup:
	@echo "Copying files from shared and hbook"
	@./shared/scripts/copy_files.sh shared
ifdef G4_USE_HBOOK
	@./shared/scripts/copy_files.sh hbook
endif
        
clean_setup:
	@echo "Removing files copied from shared and hbook"
	@./shared/scripts/clean_files.sh shared
ifdef G4_USE_HBOOK
	@./shared/scripts/clean_files.sh hbook
endif

include hbook/GNUmakefile.tools_hbook

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += -I./shared/include

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

histclean:
	rm ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/HistoManager.o
	rm ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/ExG4HbookAnalysisManager.o
