# $Id: GNUmakefile,v 1.1 2004-05-25 08:52:47 guatelli Exp $
# --------------------------------------------------------------

name := Brachy
G4TARGET := $(name)
G4EXLIB := true
G4LIB_BUILD_SHARED = 1

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
ifdef G4ANALYSIS_USE  
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib` 
endif


