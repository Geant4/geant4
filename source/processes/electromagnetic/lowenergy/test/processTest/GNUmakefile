# $Id: GNUmakefile,v 1.4 2001-10-29 12:03:15 pia Exp $
# --------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif
name := $(G4TARGET)
G4EXLIB := true
GENERAL=1

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

.PHONY: all
all: lib bin

#include $(G4INSTALL)/config/architecture.gmk

# NEW: this requires the env-var LHCXX_REL_DIR to be set up correctly,
# see the first lines in G4ComptonTest.cc
#
# C_FLAGS += -I$(LHCXX_REL_DIR)/include
# LDFLAGS  += -L$(LHCXX_REL_DIR)/lib
# general stuff:
# LOADLIBS += -lHepUtilities -lAIDA_Annotation
# add for Histograms:
# LOADLIBS += -lAIDA_HBook 
# add for Ntuples:
# LOADLIBS += -lAIDA_Ntuple -lNtupleTag -lNtupleTag_HBook 
# ... and finally system stuff and to deal with the "Fortran legacy":
# LOADLIBS += -lCHBook -ldl -lg2c 

include $(G4INSTALL)/config/binmake.gmk
