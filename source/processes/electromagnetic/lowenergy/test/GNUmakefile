# $Id: GNUmakefile,v 1.16 2003-11-10 22:31:33 silvarod Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------


ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all
all: bin

#G4EXEC_BUILD = true

ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif


include $(G4INSTALL)/config/architecture.gmk

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
endif

include $(G4INSTALL)/config/binmake.gmk




