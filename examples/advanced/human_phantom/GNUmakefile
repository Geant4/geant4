# --------------------------------------------------------------
# GNUmakefile for examples module.
# Author: Giorgio Guerrieri
# Date:  14/06/2005.
# --------------------------------------------------------------

name := phantom
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ifdef G4ANALYSIS_USE  
  CPPFLAGS  += `aida-config --include`
  EXTRALIBS += `aida-config --lib`
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
