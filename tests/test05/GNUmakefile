# $Id: GNUmakefile,v 1.2 1999-12-08 16:30:50 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test05
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all tools
 
tools:
	cd ../tools/lib; $(MAKE)
	$(MAKE) all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk
