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
ifdef ULTRA_MIRROR_USE
  CPPFLAGS += -DULTRA_MIRROR_USE
endif
ifdef ULTRA_GROUND_USE
  CPPFLAGS += -DULTRA_GROUND_USE
endif


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk


