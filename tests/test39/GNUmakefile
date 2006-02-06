# $Id: GNUmakefile,v 1.1 2006-02-06 09:16:43 mkossov Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test39
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../
endif

.PHONY: all tools
 
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

# ---Temporary--- external CHIPS processes
###CPPFLAGS +=-I$(G4BASE)/processes/photolepton_hadron/cross-sections/include \
###           -I$(G4BASE)/processes/photolepton_hadron/processes/include
###EXTRALIBS += -lG4photoleptohadr_cs -lG4photoleptohadr_proc


include $(G4INSTALL)/config/binmake.gmk

