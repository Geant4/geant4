# $Id: GNUmakefile,v 1.5 2002-09-17 13:59:14 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exampleB01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif


.PHONY: all
all: lib  bin


include $(G4INSTALL)/config/architecture.gmk


CPPFLAGS += -I$(HOME)/G4scoring/scoring/include



LDFLAGS += -L$(G4WORKDIR)/tmp/$(G4SYSTEM)/scoring 





include $(G4INSTALL)/config/binmake.gmk
MYLIBS :=  -lexampleB01 -lscoring $(LDLIBS)
LDLIBS = $(MYLIBS)

