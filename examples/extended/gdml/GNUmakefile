# $Id: GNUmakefile,v 1.6 2006/07/24 15:32:13 witoldp Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# Author: Radovan Chytracek
# Date:   04/12/2004.
# --------------------------------------------------------------

name := g4_gdml_read_write
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

CPPFLAGS    += \
		-I$(GDML_INC) \
	       	-I$(XERCESCROOT)/include

EXTRALIBS   += \
               -L$(GDML_LIB) \
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
ifndef XERCESCROOT	
	@echo XERCESCROOT not defined!
endif
ifndef GDML_INC
	@echo GDML_INC not defined!
endif
ifndef GDML_LIB
	@echo GDML_LIB not defined!
endif


include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
