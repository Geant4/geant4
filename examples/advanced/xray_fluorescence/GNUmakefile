# $Id: GNUmakefile,v 1.14 2003-07-17 14:16:05 mantero Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

# G4LIB_BUILD_SHARED = 1

ifndef G4TARGET
name := XrayFluo
G4TARGET := $(name)
endif

G4EXLIB := true

G4ANALYSIS_USE = 0

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all visclean
all:  lib bin 
include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib`
endif

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

