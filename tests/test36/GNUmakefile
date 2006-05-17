# $Id: GNUmakefile,v 1.1 2006-05-17 16:09:24 grichine Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test36
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../
endif

.PHONY: all
all: lib bin

ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif

include $(G4INSTALL)/config/architecture.gmk

# for the aida-config command see the README file
ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
endif

include $(G4INSTALL)/config/binmake.gmk

histclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Em10RunAction.o

visclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/Em10VisManager.o

