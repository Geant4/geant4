# $Id: GNUmakefile,v 1.1 1999-01-08 16:31:50 gunter Exp $
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
INCFLAGS := -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/graphics_reps/include

LDLIBS   := -lG4volumes -lG4csg -lG4geometrymng -lG4volumes -lG4geometrymng \
            -lG4global \
            -lG4graphics_reps

include $(G4INSTALL)/config/binmake.gmk
