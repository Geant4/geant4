# $Id: GNUmakefile,v 1.1 2002-11-18 10:50:55 jwellisc Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test25
G4TARGET := $(name)
G4EXLIB := true
EXTRALIBS :=-pg 

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all tools
 
tools:
	cd ../tools/lib; $(MAKE)
	$(MAKE) all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk
