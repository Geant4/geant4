# $Id: GNUmakefile,v 1.4 2001-09-18 12:34:49 pfeiffer Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../..
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
LOADLIBS += -l AIDA_Ntuple -lNtupleTag -lNtupleTag_HBook 
# ... and finally system stuff and to deal with the "Fortran legacy":
LOADLIBS += -lCHBook -ldl -lg2c 

# old:
#LDFLAGS  += -L/cern/pro/lib
#LOADLIBS += -lpacklib $(FCLIBS)

include $(G4INSTALL)/config/binmake.gmk


