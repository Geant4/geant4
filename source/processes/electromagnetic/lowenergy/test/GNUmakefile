# $Id: GNUmakefile,v 1.3 2001-09-14 16:52:20 pia Exp $
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

C_FLAGS += -I$(LHCXX_REL_DIR)/include
LDFLAGS  += -L$(LHCXX_REL_DIR)/lib
LOADLIBS += -lHepUtilities -lAIDA_HBook -lCHBook -lAIDA_Annotation -ldl -lg2c 

#LDFLAGS  += -L/cern/pro/lib
#LOADLIBS += -lpacklib $(FCLIBS)

include $(G4INSTALL)/config/binmake.gmk



