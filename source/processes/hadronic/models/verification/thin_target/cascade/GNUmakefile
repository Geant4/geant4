# $Id: GNUmakefile,v 1.2 2003-11-14 17:02:42 vnivanch Exp $
# --------------------------------------------------------------

name := cascade
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../../../../
endif

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


