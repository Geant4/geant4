# $Id: GNUmakefile,v 1.17 2006-03-06 16:41:11 silvarod Exp $
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

G4ANALYSIS_USE = 1

ifdef G4ANALYSIS_USE
  CPPFLAGS += -DG4ANALYSIS_USE
endif


include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
  LOADLIBS += `aida-config --lib`
#   CPPFLAGS += -I/afs/cern.ch/sw/geant4/dev/PI/1.3.3-lite/redhat73_gcc323/include
#   LDFLAGS  += -L/afs/cern.ch/sw/geant4/dev/PI/1.3.3-lite/redhat73_gcc323/lib
endif







