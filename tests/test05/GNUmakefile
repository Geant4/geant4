# $Id: GNUmakefile,v 1.3 2001-11-29 16:50:58 gcosmo Exp $
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

.PHONY: all tools
 
tools:
	cd ../tools/lib; $(MAKE)
	$(MAKE) all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk
