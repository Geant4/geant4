# $Id: GNUmakefile,v 1.2 2006-11-16 15:02:21 vnivanch Exp $
# --------------------------------------------------------------

name := test35
G4TARGET := test35
G4EXLIB  := true


.PHONY: all tools


all: lib bin


ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LOADLIBS += `aida-config --lib`
endif


include $(G4INSTALL)/config/binmake.gmk


