# $Id: GNUmakefile,v 1.12 2002-06-16 14:57:45 mantero Exp $
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

CPPFLAGS += -I${LHCXX_REL_DIR}/include -I/afs/cern.ch/sw/contrib/AIDA/2.2/cpp
LOADLIBS  += -L/usr/local/gcc-alt-2.95.2



visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

