# $Id: GNUmakefile,v 1.5 2002-11-08 14:47:31 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleB02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin


CPPFLAGS +=

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += `aida-config -i`
LDFLAGS += `aida-config -l`
