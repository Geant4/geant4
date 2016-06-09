# $Id: GNUmakefile,v 1.1 2008/11/03 11:48:35 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := pythia6_decayer
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all

all: pythia6 lib bin

EXTRALIBS = $(G4TMPDIR)/libPythia6.so

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

