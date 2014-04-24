# $Id$
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../../../
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk
	



ifdef G4ANALYSIS_USE_ROOT
  CPPFLAGS += `root-config --cflags --ldflags --libs`
  LOADLIBS += `root-config --ldflags --libs`
endif

include $(G4INSTALL)/config/binmake.gmk
    










