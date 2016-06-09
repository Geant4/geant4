# $Id: GNUmakefile,v 1.17 2004/05/29 10:02:08 mantero Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

# G4LIB_BUILD_SHARED = 1

ifndef G4TARGET
name := XrayFluo
G4TARGET := $(name)
endif

G4EXLIB := true

G4ANALYSIS_USE = 1

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all visclean
all:  bin #lib 
include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
# LOADLIBS += `aida-config --lib`
endif


# As the dependency is checked before "all" is reached, we need
# to do some tricks here ...

# to make simulation working uncomment the line below
#CPPFLAGS += -I${LHCXX_REL_DIR}/include -I/afs/cern.ch/sw/contrib/AIDA/2.2/cpp
#LOADLIBS  += -L/users/mantero/lib -L/users/mantero/lib/gcc-lib/2.96 -L/users/mantero/lib/gcc-lib/i686-pc-linux-gnu/2.95.2
LOADLIBS  += -L/users/mantero/lib -lf2c -L/usr/local/lib/ -lg2c
#LOADLIBS  += -L/users/mantero/lib/gcc-lib/i686-pc-linux-gnu/2.95.2 -lg2c

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

