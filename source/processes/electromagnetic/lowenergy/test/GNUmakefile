# $Id: GNUmakefile,v 1.6 2001-11-06 09:46:20 pia Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif
name:= $(G4TARGET)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

#
# NEW: this requires the env-var LHCXX_REL_DIR to be set up correctly,
# see the first lines in G4ComptonTest.cc
#
CPPFLAGS += -I$(LHCXX_REL_DIR)/include
LDFLAGS  += -L$(LHCXX_REL_DIR)/lib
# general stuff:
LOADLIBS += -lHepUtilities -lAIDA_Annotation
# add for Histograms:
LOADLIBS += -lAIDA_HBook 
# add for Ntuples:
LOADLIBS += -lAIDA_Ntuple -lNtupleTag -lNtupleTag_HBook 
# ... and finally system stuff and to deal with the "Fortran legacy":
LOADLIBS += -lCHBook -ldl -lg2c 

# old:
#LDFLAGS  += -L/cern/pro/lib
#LOADLIBS += -lpacklib $(FCLIBS)

include $(G4INSTALL)/config/binmake.gmk


