# $Id: GNUmakefile,v 1.2 2004-06-01 07:05:03 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for LXe example.
# --------------------------------------------------------------

name := LXe
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk
