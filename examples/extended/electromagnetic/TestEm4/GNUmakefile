# $Id: GNUmakefile,v 1.6 2002-05-29 15:32:26 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm4
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

###G4NOHIST := true
ifdef G4NOHIST
  CPPFLAGS += -DG4NOHIST
else
  # for the aida-config command see the README file
  CPPFLAGS += `aida-config --useHBook --include`
  LDFLAGS  += `aida-config --useHBook --lib`
endif

include $(G4INSTALL)/config/binmake.gmk

