# $Id: GNUmakefile,v 1.15 2003/12/04 11:01:47 guatelli Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayFluo
G4TARGET := $(name)
G4EXLIB := true

G4ANALYSIS_USE = 1

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

.PHONY: all visclean
all:  lib bin 

include $(G4INSTALL)/config/binmake.gmk

# Copy the XrayFluoAnalysisManage files according to which version
# is selected.
# As the dependency is checked before "all" is reached, we need
# to do some tricks here ...

ifdef G4ANALYSIS_USE  
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib` 
endif

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

