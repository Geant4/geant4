# $Id: GNUmakefile,v 1.2 2002-12-09 11:25:10 gcosmo Exp $
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

.PHONY: all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk
