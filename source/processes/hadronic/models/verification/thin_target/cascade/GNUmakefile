# $Id: GNUmakefile,v 1.5 2006-03-03 14:26:47 vnivanch Exp $
# --------------------------------------------------------------

G4TARGET := cascade_test
name := G4TARGET
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


