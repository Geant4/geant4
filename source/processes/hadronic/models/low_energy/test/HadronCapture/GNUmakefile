# $Id: GNUmakefile,v 1.1 2008-08-27 22:12:44 dennis Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := hadroncapture
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

# Root (exlude libNew and libpthread from library list)
ROOTINC       = -I$(ROOTSYS)/include
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit -lHtml
ROOTLIBS      := $(filter-out -lNew,$(ROOTLIBS))
ROOTLIBS      := $(filter-out -lpthread,$(ROOTLIBS))
  
# Extra flags for G4
  
CPPFLAGS += $(ROOTINC)
  
LDLIBS   += $(ROOTLIBS)
#
CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

