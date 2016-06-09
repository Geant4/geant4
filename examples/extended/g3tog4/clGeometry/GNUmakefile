# $Id: GNUmakefile,v 1.2 2006-03-14 10:13:55 vnivanch Exp $
# ----------------------------------------------------------------

name := clGeometry
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk
