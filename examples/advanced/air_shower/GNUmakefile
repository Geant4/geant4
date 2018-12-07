# --------------------------------------------------------------
# GNUmakefile for ULTRA-based air shower Advanced Example.  
# --------------------------------------------------------------

name := Ultra
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef SPACE_VERBOSE
  CPPFLAGS += -DULTRA_VERBOSE_USE
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk


