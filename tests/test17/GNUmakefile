# $Id: GNUmakefile,v 1.4 2000-05-26 15:50:37 pia Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test17
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../
endif

.PHONY: all tools
 
tools:
	cd ../tools/lib; $(MAKE)
	$(MAKE) all

all: lib bin

include $(G4INSTALL)/config/architecture.gmk


# ifeq ($(G4SYSTEM),Linux-g++)
# LDFLAGS  += -L/afs/cern.ch/asis/i386_redhat51/cern/pro/lib
# LOADLIBS += -lpacklib $(FCLIBS)
# else
# LDFLAGS  += -L/cern/pro/lib
# LOADLIBS += -lpacklib $(FCLIBS)
# endif

include $(G4INSTALL)/config/binmake.gmk


