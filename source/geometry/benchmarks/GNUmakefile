# $Id: GNUmakefile,v 1.4 2002-01-09 16:17:55 gcosmo Exp $
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
INCFLAGS := -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/solids/verification/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/magneticfield/include \
	    -I$(G4BASE)/materials/include \
            -I$(G4BASE)/graphics_reps/include \
	    -I$(G4BASE)/intercoms/include

LDLIBS   := \
	    -lG4csg \
	    -lG4volumes \
	    -lG4magneticfield \
            -lG4geomver \
	    -lG4geometrymng \
	    -lG4materials \
            -lG4graphics_reps \
	    -lG4intercoms \
            -lG4globman
ifdef G4GEOMETRY_VERBOSE
  CPPFLAGS += -DG4GEOMETRY_VERBOSE
endif

include $(G4INSTALL)/config/binmake.gmk
