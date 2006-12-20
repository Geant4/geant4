# $Id: GNUmakefile,v 1.4 2006-12-20 13:50:24 flei Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := exGPS
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

#G4VIS_BUILD_OPENGLX_DRIVER = 0
#G4VIS_USE_OPENGLX = 0
#LDFLAGS += -static

ifdef G4ANALYSIS_USE_ROOT
CPPFLAGS += -DG4ANALYSIS_USE_ROOT -DG4ANALYSIS_USE
CPPFLAGS += -D_REENTRANT -I$(ROOTSYS)/include

EXTRALIBS += -L$(ROOTSYS)/lib \
-lCore -lCint \
-lHist -lGraf -lGraf3d -lGpad \
-lTree -lRint -lPostscript \
-lMatrix -lPhysics \
-lm -ldl -lpthread -rdynamic

endif

ifdef G4ANALYSIS_USE_AIDA
   CPPFLAGS += -DG4ANALYSIS_USE_AIDA -DG4ANALYSIS_USE
   CPPFLAGS += `aida-config --include`
   LDLIBS += `aida-config --lib`
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

#LDFLAGS += -static
