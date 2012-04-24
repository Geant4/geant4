# $Id:$
# --------------------------------------------------------------
# GNUmakefile for beamTest.
# --------------------------------------------------------------

name := test73
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all:  lib bin

ifndef G4ANALYSIS_USEROOT
   $(warning This test need ROOT to run properly, code will not run, \
     set G4ANALYSIS_USEROOT environment variable to 1)
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USEROOT
   EXTRALIBS += `root-config --libs`
   CPPFLAGS+=  `root-config --cflags`
   CPPFLAGS+= -DG4ANALYSIS_USEROOT
endif

include $(G4INSTALL)/config/binmake.gmk

