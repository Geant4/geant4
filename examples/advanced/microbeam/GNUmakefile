# -------------------------------------------------------------------
# $Id: GNUmakefile,v 1.3 2006-04-10 14:47:31 sincerti Exp $
# -----------
# GNUmakefile 
# -----------

name := Microbeam

G4DEBUG := 0
G4TARGET := $(name)
G4EXLIB := true
G4VERBOSE :=1

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
