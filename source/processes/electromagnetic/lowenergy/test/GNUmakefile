# $Id: GNUmakefile,v 1.19 2009-05-16 11:07:23 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

G4EXEC_BUILD := true

ifndef G4INSTALL
  G4INSTALL = ../../../../../
endif


ifdef G4ANALYSIS_USE_ROOT
  CPPFLAGS += -DG4ANALYSIS_USE_ROOT
endif

ifdef X11LIBS 
  LDFLAGS += $(X11LIBS)
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE_ROOT
  CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
  LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)
endif

include $(G4INSTALL)/config/binmake.gmk

#
