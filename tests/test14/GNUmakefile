# $Id: GNUmakefile,v 1.3 1999-12-08 16:30:55 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test14
G4TARGET := $(name)
G4EXLIB := true
GENERAL=1

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all tools
 
tools:
	cd ../tools/lib; $(MAKE)
	$(MAKE) all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk
