# $Id: GNUmakefile,v 1.4 2003/05/28 08:12:48 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := cadex
G4TARGET := $(name)
G4EXLIB := true
G4USE_STEP := 1

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

# To avoid incompatibilities with ANSI standard treatment of strings, the
# NIST STEP reader on gcc-2.95.X and later must be compiled in non-ANSI mode.
# NOTE: correct behavior of the STEP reader and interface are guaranteed ONLY
#       if using egcs-1.1.2 compiler or old non-ANSI compilers.
#
ifneq (,$(findstring g++,$(G4SYSTEM)))
  ansi_options = -Wall -ansi -pedantic
  NCXX = $(filter-out $(ansi_options),$(CXXFLAGS))
  NCPP = $(CPPFLAGS)
  CXXFLAGS := $(NCXX)
  CPPFLAGS := $(NCPP) -w
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
