# $Id: GNUmakefile,v 1.1 1999-01-08 16:32:10 gunter Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables for binmake.gmk.
#
INCFLAGS := -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/materials/include

LDLIBS   := \
            -lG4materials \
            -lG4processes \
            -lG4track \
            -lG4particles \
            -lG4intercoms \
            -lG4global \
            -lG4intercoms

include $(G4INSTALL)/config/binmake.gmk
