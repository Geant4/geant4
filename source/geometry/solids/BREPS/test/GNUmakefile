# $Id: GNUmakefile,v 1.2 1999-06-07 17:44:28 japost Exp $
# ------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables from GNUmakefile.architecture.
INCFLAGS := \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
	    -I$(G4BASE)/geometry/solids/STEP/include \
	    -I$(G4BASE)/geometry/solids/STEPinterface/include \
	    -I$(G4BASE)/geometry/solids/BREPS/include \
	    -I$(G4BASE)/materials/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include
# Try to move STEP interface tests to that test dir --> remove dir here
# Try to delete materials from test -> remove here

LDLIBS   := \
	    -lG4stepinterface \
	    -lG4brep \
	    -lG4step \
	    -lG4csg \
	    -lG4geometrymng \
            -lG4graphics_reps \
            -lG4globman

include $(G4INSTALL)/config/binmake.gmk
