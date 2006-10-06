# $Id: GNUmakefile,v 1.6 2006-10-06 13:10:30 gunter Exp $
# --------------------------------------------------------------

G4TARGET := cascade_test
name := $(G4TARGET)
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


