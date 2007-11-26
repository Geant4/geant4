# $Id: GNUmakefile,v 1.1 2007-11-26 16:22:56 grichine Exp $
# --------------------------------------------------------------
# GNUmakefile for LXe example.
# --------------------------------------------------------------

name := test42
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk
