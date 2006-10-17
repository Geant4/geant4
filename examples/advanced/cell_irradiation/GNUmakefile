# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.1 2006-10-17 17:57:24 mascial Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := cell
G4TARGET := $(name)
G4EXLIB := true

G4ANALYSIS_USE = 1

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE  
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib` 
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk


