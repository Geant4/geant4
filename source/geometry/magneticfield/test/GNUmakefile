# $Id: GNUmakefile,v 1.2 1999-07-06 20:25:41 japost Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable TESTTARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables from binmake.gmk.
#
INCFLAGS := -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/graphics_reps/include

OLD_LDLIBS   := -lG4geometry  \
            -lG4graphics_reps \
            -lG4global

LDLIBS   := -lG4csg -lG4volumes -lG4geometrymng -lG4magneticfield \
            -lG4graphics_reps \
            -lG4globman 

include $(G4INSTALL)/config/binmake.gmk
