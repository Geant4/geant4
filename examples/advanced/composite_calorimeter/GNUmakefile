# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := CompositeCalorimeter
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef debug
  CPPFLAGS += -Ddebug
endif

.PHONY: all
 
all: lib bin

include $(G4INSTALL)/config/binmake.gmk



