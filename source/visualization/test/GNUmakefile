# $Id: GNUmakefile,v 1.3 2001-08-17 22:56:39 johna Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test19
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

CPPFLAGS += -I$(G4BASE)/visualization/XXX/include
EXTRALIBS += -lG4visXXX

include $(G4INSTALL)/config/binmake.gmk
