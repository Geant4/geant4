# $Id: GNUmakefile,v 1.2 2001-10-15 17:36:13 pia Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# G4NOHIST variable is checked.  Vladimir Ivanchenko, 2000/01/18
# --------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

#
# NEW: this requires the env-var LHCXX_REL_DIR to be set up correctly,
# see the first lines in G4ComptonTest.cc
#
C_FLAGS += -I$(LHCXX_REL_DIR)/include
LDFLAGS  += -L$(LHCXX_REL_DIR)/lib
# general stuff:
LOADLIBS += -lHepUtilities -lAIDA_Annotation
# add for Histograms:
LOADLIBS += -lAIDA_HBook 
# add for Ntuples:
LOADLIBS += -lAIDA_Ntuple -lNtupleTag -lNtupleTag_HBook 
# ... and finally system stuff and to deal with the "Fortran legacy":
LOADLIBS += -lCHBook -ldl -lg2c 

include $(G4INSTALL)/config/binmake.gmk




