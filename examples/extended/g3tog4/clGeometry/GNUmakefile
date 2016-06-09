# $Id: GNUmakefile,v 1.2 2006-03-14 10:13:55 vnivanch Exp $
# ----------------------------------------------------------------

name := clGeometry
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: setup clean_setup all
all: lib bin

setup:
	@echo "Copying files from common"
	@$(G4INSTALL)/examples/extended/common/scripts/copy_files.sh

clean_setup:
	@echo "Removing files copied from common"
	@$(G4INSTALL)/examples/extended/common/scripts/clean_files.sh

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk
