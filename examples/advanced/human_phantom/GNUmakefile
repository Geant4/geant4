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
   GDML_BASE := /local/GDML_2_8_0
endif

CPPFLAGS    += \
               -I$(GDML_BASE)/include \
               -I$(XERCESCROOT)/include 

ifdef G4ANALYSIS_USE  
CPPFLAGS    += \
	         `aida-config --include`
endif

EXTRALIBS   += \
               -L$(GDML_BASE)/$(SCRAM_ARCH)/lib \
               -lG4Processor \
               -lG4Subscribers \
               -lG4Processor \
	       -lG4Evaluator\
               -lG4Writer \
               -L$(GDML_BASE)/$(SCRAM_ARCH)/lib \
               -lSaxana \
               -lProcesses \
               -lWriter \
               -L$(XERCESCROOT)/lib -lxerces-c 
           
ifdef G4ANALYSIS_USE  
EXTRALIBS    += \
	         `aida-config --lib`
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
