# $Id: GNUmakefile,v 1.4 2002-12-09 11:25:09 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test05
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

# To avoid incompatibilities with ANSI standard treatment of strings, the
# NIST STEP reader on gcc-2.95.X and later must be compiled in non-ANSI mode.
#
ifneq (,$(findstring Linux-g++,$(G4SYSTEM)))
  ansi_options = -Wall -ansi -pedantic
  NCXX = $(filter-out $(ansi_options),$(CXXFLAGS))
  NCPP = $(CPPFLAGS)
  CXXFLAGS := $(NCXX)
  CPPFLAGS := $(NCPP) -w
endif

.PHONY: all
 
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
