# $Id: GNUmakefile,v 1.4 2006/05/30 07:22:42 mkossov Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := PhotIn
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

# New Lists definition is necessary only for comparison with the new lists
###include hadronic_lists.gmk
# New Lists definition is necessary only while it is out of the G4 kernel
###include new_phys_lists.gmk

#G4ANALYSIS_USE:= 1

.PHONY: all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk

##ifdef G4ANALYSIS_USE
##  CPPFLAGS += `aida-config --include`
##  LDFLAGS  += `aida-config --lib`
##  LOADLIBS += `aida-config --lib`
##endif

#visclean:
#	rm -f g4*.prim g4*.eps g4*.wrl
#	rm -f .DAWN_*

