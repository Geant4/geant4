# $Id: GNUmakefile,v 1.2 2003-04-01 10:59:39 gcosmo Exp $
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
LDLIBS   := -lG4globman

include $(G4INSTALL)/config/binmake.gmk
