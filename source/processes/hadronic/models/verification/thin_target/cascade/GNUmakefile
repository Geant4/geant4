# $Id: GNUmakefile,v 1.3 2003-11-14 17:39:20 vnivanch Exp $
# --------------------------------------------------------------

name := cascade
G4TARGET := cascade
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


