# -------------------------------------------------------------------
# $Id: GNUmakefile,v 1.4.6.1 2010/03/18 11:14:26 gcosmo Exp $
# -----------
# GNUmakefile 
# -----------

name := Microbeam

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
