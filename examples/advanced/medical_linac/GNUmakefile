# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.2 2005/07/03 23:27:36 mpiergen Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := MedLinac

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin
include $(G4INSTALL)/config/binmake.gmk


ifdef G4ANALYSIS_USE  
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib` 
endif

ifdef G4STORE_DATA
  CPPFLAGS += -DG4STORE_DATA
endif

