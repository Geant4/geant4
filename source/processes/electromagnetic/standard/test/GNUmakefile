# $Id: GNUmakefile,v 1.5 2006-10-20 08:59:50 vnivanch Exp $
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

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)

include $(G4INSTALL)/config/binmake.gmk
