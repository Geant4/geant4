# $Id: GNUmakefile,v 1.4 2006-06-01 10:15:14 witoldp Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# Author: Radovan Chytracek
# Date:   04/12/2004.
# --------------------------------------------------------------

name := g4_gdml_read_write
G4TARGET := $(name)
G4EXLIB := true

ifndef SCRAM_ARCH
  SCRAM_ARCH := slc3_ia32_gcc323
endif

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ifndef GDML_BASE
  GDML_BASE := /afs/cern.ch/sw/lcg/app/releases/GDML/GDML_2_6_0
endif

ifndef XERCESCROOT
  XERCESCROOT := /afs/cern.ch/sw/lcg/external/XercesC/2.7.0/${SCRAM_ARCH}
endif

CPPFLAGS    += \
		-I$(GDML_BASE)/include \
	       	-I$(XERCESCROOT)/include

EXTRALIBS   += \
               -L$(GDML_BASE)/$(SCRAM_ARCH)/lib \
               -lG4Evaluator \
               -lG4Processor \
               -lG4Subscribers \
               -lG4Writer \
               -lSaxana \
               -lProcesses \
               -lWriter \
               -L$(XERCESCROOT)/lib -lxerces-c 

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
