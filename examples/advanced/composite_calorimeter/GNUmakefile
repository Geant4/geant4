# $Id: GNUmakefile,v 1.9 2003/12/04 11:20:34 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := CompositeCalorimeter
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../
endif

include $(G4INSTALL)/config/architecture.gmk

include hadronic_lists.gmk

.PHONY: all
 
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE  
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
  LOADLIBS += `aida-config --lib`
endif

