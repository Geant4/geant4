# $Id: GNUmakefile,v 1.4 2004-02-10 18:07:27 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := CrossSectionTest
endif

G4EXEC_BUILD := true


ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
endif

include $(G4INSTALL)/config/binmake.gmk


