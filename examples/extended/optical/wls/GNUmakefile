# $Id: GNUmakefile,v 1.1 2009/11/25 23:36:12 gum Exp $
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

ifdef MAKEHBOOK
 CPPFLAGS+=-I$(G4INSTALL)/tests/tools/include -DMAKEHBOOK -DMAKEBOOK
 EXTRALIBS+=-L$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM) -lG4TestTool
 EXTRALIBS+=-L/cern/pro/lib -lpacklib -lnsl
 ifeq ($(G4SYSTEM),Linux-g++)
   EXTRALIBS+=-lg2c -lm -lgcc -lc -lgcc
 endif
endif

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk
