# $Id$
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
            -I$(G4BASE)/global/HEPRandom/include
LDLIBS   := -lG4heprandom -lG4globman
ifdef G4LIB_USE_CLHEP
  LDLIBS += -lG4clhep
endif

include $(G4INSTALL)/config/binmake.gmk
