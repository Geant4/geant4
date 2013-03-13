# $Id$
# --------------------------------------------------------------
# GNUmakefile for Bertini test executables
# --------------------------------------------------------------

name := test75
G4TARGET := $(name)

#include $(G4INSTALL)/config/binmake.gmk
include ../test48/binmake.gmk

# Root (exlude libNew and libpthread from library list)
ROOTINC  := -I$(ROOTSYS)/include
ROOTLIBS := $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit -lHtml
ROOTLIBS := $(filter-out -lNew,$(ROOTLIBS))
ROOTLIBS := $(filter-out -lpthread,$(ROOTLIBS))

# Extra flags for G4
CPPFLAGS += $(ROOTINC)
LDLIBS   += $(ROOTLIBS)
