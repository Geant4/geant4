# $Id: GNUmakefile,v 1.1 2006-01-11 15:43:52 mkossov Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := AnaEx02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

#ifdef G4ANALYSIS_USE
# CPPFLAGS += `aida-config --include`
# LDFLAGS  += `aida-config --lib`
# LOADLIBS += `aida-config --lib`
#endif

