# $Id: GNUmakefile 72244 2013-07-12 08:49:56Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := pythia6_decayer
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: setup clean_setup all

all: pythia6 lib bin

EXTRALIBS = $(G4TMPDIR)/libPythia6.so

setup:
	@echo "Copying files from common"
	@$(G4INSTALL)/examples/extended/common/scripts/copy_files.sh

clean_setup:
	@echo "Removing files copied from common"
	@$(G4INSTALL)/examples/extended/common/scripts/clean_files.sh

include $(G4INSTALL)/config/binmake.gmk

CCFLAGS += -c
FCFLAGS += -c

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

pythia6: $(G4TMPDIR)/libPythia6.so

$(G4TMPDIR)/libPythia6.so: $(G4TMPDIR)/pythia6.o $(G4TMPDIR)/pythia6_common_address.o
	$(FC) -shared -Wl,-soname,libPythia6.so -o $(G4TMPDIR)/libPythia6.so  $(G4TMPDIR)/pythia6.o $(G4TMPDIR)/pythia6_common_address.o 
 
$(G4TMPDIR)/pythia6.o:
	$(FC) $(FCFLAGS) $(PYTHIA6)/pythia-$(PYTHIA6_VERSION).f -o $(G4TMPDIR)/pythia6.o

$(G4TMPDIR)/pythia6_common_address.o:
	$(CC) $(CCFLAGS) src/pythia6_common_address.c -o $(G4TMPDIR)/pythia6_common_address.o

