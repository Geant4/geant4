# $Id: GNUmakefile,v 1.4 1999/12/01 15:46:30 morita Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------
# You must have unique G4EXAMPLE_FDID assigned by your local 
# manager of Objectivity/DB lock server.
#

name := PersistentEx01
G4ODBMS := true
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all doall

all: cleandb headers
	@(G4ODBMS_DEPEND=1; export G4ODBMS_DEPEND; $(MAKE) doall)

doall: lib bin

include $(G4INSTALL)/config/architecture.gmk

G4OODDLX_BOOT := $(G4EXAMPLE_BOOT)

include $(G4INSTALL)/config/binmake.gmk

  CPPFLAGS += \
              -Iinclude \
              -I$(G4TMPDIR)

