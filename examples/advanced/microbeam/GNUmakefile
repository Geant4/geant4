# -------------------------------------------------------------------
# $Id: GNUmakefile,v 1.1 2006-04-04 13:02:15 sincerti Exp $
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
