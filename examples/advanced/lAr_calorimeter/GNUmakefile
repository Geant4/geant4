# $Id: GNUmakefile,v 1.5 2003/05/26 12:51:17 pmendez Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := lArCal
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

#ifdef G4ANALYSIS_USE
#  CPPFLAGS += -DG4ANALYSIS_USE
#endif


include $(G4INSTALL)/config/architecture.gmk

G4ANALYSIS_USE:= true
#ifdef G4ANALYSIS_USE
#  CPPFLAGS += -DG4ANALYSIS_USE
#endif

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS += `aida-config --lib`
endif


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*












