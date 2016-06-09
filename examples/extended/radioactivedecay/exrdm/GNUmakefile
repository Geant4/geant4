# $Id: GNUmakefile,v 1.4 2006/12/21 17:52:39 flei Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := exRDM
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

ifdef G4ANALYSIS_USE_ROOT
CPPFLAGS += -DG4ANALYSIS_USE_ROOT 
CPPFLAGS += -D_REENTRANT -I$(ROOTSYS)/include

EXTRALIBS += -L$(ROOTSYS)/lib \
-lCore -lCint \
-lHist -lGraf -lGraf3d -lGpad \
-lTree -lRint -lPostscript \
-lMatrix -lPhysics \
-lm -ldl -lpthread -rdynamic

endif

ifdef G4ANALYSIS_USE_AIDA
   CPPFLAGS += -DG4ANALYSIS_USE_AIDA
   CPPFLAGS += `aida-config --include`
   LDLIBS += `aida-config --lib`
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

