# $Id: GNUmakefile,v 1.7 2007-11-16 15:21:38 ztorzsok Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.
# Author: Radovan Chytracek
# Date:   04/12/2004.
# --------------------------------------------------------------

name := load_gdml
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

CPPFLAGS += -I$(XERCESCROOT)/include \
            -DG4LIB_USE_XERCESC

EXTRALIBS += -lG4gdml \
             -L$(XERCESCROOT)/lib -lxerces-c 

.PHONY: all
all: lib bin
ifndef XERCESCROOT	
	@echo XERCESCROOT not defined!
endif

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
