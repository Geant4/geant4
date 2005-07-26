# $Id: GNUmakefile,v 1.1 2005-07-26 13:41:22 gguerrie Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# Author: Giorgio Guerrieri
# Date:  14/06/2005.
# --------------------------------------------------------------

name := phantom
G4TARGET := $(name)
G4EXLIB := true

ifndef SCRAM_ARCH
   SCRAM_ARCH := slc3_ia32_gcc323
endif

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ifndef GDML_BASE
#  GDML_BASE := /afs/cern.ch/sw/lcg/app/releases/GDML/pro
   GDML_BASE := /home3/gguerrie/devel/GDML_2_1_0
endif

CPPFLAGS    += \
               -I$(GDML_BASE)/include \
               -I$(GDML_BASE)/geant4/include \
               -I$(XERCESCROOT)/include \
               `aida-config --include`

EXTRALIBS   += \
               -L$(GDML_BASE)/geant4/$(SCRAM_ARCH)/lib \
               -llcg_G4Processor \
               -llcg_G4Subscribers \
               -llcg_G4Processor \
               -llcg_G4Writer \
               -L$(GDML_BASE)/$(SCRAM_ARCH)/lib \
               -llcg_Saxana \
               -llcg_Processes \
               -llcg_Writer \
               -L$(XERCESCROOT)/lib -lxerces-c \
               `aida-config --lib`

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
