# $Id: GNUmakefile,v 1.1 1999-01-08 16:32:06 gunter Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

# Override some variables for binmake.gmk.
#
# INCFLAGS :=
# LDLIBS   :=

include $(G4INSTALL)/config/binmake.gmk
