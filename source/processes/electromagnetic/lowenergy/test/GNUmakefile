# $Id: GNUmakefile,v 1.1 1999-06-24 19:30:17 aforti Exp $
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

# Override some variables for binmake.gmk.
#
include $(G4INSTALL)/config/binmake.gmk
LDLIBS += `cernlib packlib` $(FCLIBS)



