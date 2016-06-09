# $Id: GNUmakefile,v 1.4 2010-11-08 10:38:44 maire Exp $
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
all: hbook lib bin

setup:
	@echo "Copying files from common"
	@$(G4INSTALL)/examples/extended/common/scripts/copy_files.sh

clean_setup:
	@echo "Removing files copied from common"
	@$(G4INSTALL)/examples/extended/common/scripts/clean_files.sh

# HBOOK support
#
#### G4_USE_HBOOK := true
include GNUmakefile.tools_hbook

include $(G4INSTALL)/config/binmake.gmk


visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

histclean:
	rm ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/HistoManager.o
	rm ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/ExG4HbookAnalysisManager.o
