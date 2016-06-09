# $Id: GNUmakefile,v 1.3 2004/12/06 14:26:14 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# Author: Radovan Chytracek
# Date:   04/12/2004.
# --------------------------------------------------------------

name := g4_gdml_read_write
G4TARGET := $(name)
G4EXLIB := true

ifndef SCRAM_ARCH
  SCRAM_ARCH := rh73_gcc323
endif

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ifndef GDML_BASE
  GDML_BASE := /afs/cern.ch/sw/lcg/app/releases/GDML/pro
endif

CPPFLAGS    += \
               -I$(GDML_BASE)/include \
               -I$(GDML_BASE)/geant4/include \
               -I$(XERCESCROOT)/include

EXTRALIBS   += \
               -L$(GDML_BASE)/geant4/$(SCRAM_ARCH)/lib \
               -llcg_G4GDMLProcessor \
               -llcg_G4GDMLSubscribers \
               -llcg_G4GDMLProcessor \
               -llcg_G4GDMLWriter \
               -L$(GDML_BASE)/$(SCRAM_ARCH)/lib \
               -llcg_GDMLSaxana \
               -llcg_GDMLProcesses \
               -llcg_GDMLWriter \
               -L$(XERCESCROOT)/lib -lxerces-c

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
