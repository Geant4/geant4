# $Id: GNUmakefile,v 1.1 1999-01-07 16:08:52 gunter Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables for binmake.gmk.
#
INCFLAGS := -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPGeometry/include
LDLIBS   := -lG4global

include $(G4INSTALL)/config/binmake.gmk
