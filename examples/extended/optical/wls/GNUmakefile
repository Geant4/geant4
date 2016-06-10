# $Id: GNUmakefile 68752 2013-04-05 10:23:47Z gcosmo $
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
