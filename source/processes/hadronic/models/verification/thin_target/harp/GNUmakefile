# $Id: GNUmakefile,v 1.1 2004-03-30 10:12:25 vnivanch Exp $
# --------------------------------------------------------------

name := harp
G4TARGET := harp
G4EXLIB  := true


.PHONY: all tools


all: lib bin


ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
endif


include $(G4INSTALL)/config/binmake.gmk


