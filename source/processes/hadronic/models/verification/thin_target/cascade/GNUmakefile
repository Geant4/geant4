# $Id: GNUmakefile,v 1.4 2004-07-02 11:37:38 vnivanch Exp $
# --------------------------------------------------------------

name := cascade_test
G4TARGET := cascade_test
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


