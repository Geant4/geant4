# $Id: GNUmakefile,v 1.2 2010-01-11 16:47:50 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for WLS example.
# --------------------------------------------------------------

name := wls
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk
