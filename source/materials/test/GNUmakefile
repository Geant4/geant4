# $Id: GNUmakefile,v 1.3 2005-02-22 10:11:09 maire Exp $
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

INCFLAGS := \
	    -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/materials/include \
	    -I$(G4BASE)/intercoms/include
	    
LDLIBS   := \
            -lG4materials \
            -lG4intercoms \
            -lG4globman
	    	    
include $(G4INSTALL)/config/binmake.gmk
