# $Id: GNUmakefile,v 1.8 2002-11-08 13:58:28 ahoward Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := DMX
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

ifdef DMXENV_GPS_USE
  CPPFLAGS += -DDMXENV_GPS_USE
endif

ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif


include $(G4INSTALL)/config/architecture.gmk

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

