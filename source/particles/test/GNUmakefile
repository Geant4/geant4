# $Id: GNUmakefile,v 1.1 1999-01-08 16:32:12 gunter Exp $
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
include $(G4INSTALL)/config/G4UI_USE.gmk
include $(G4INSTALL)/config/interactivity.gmk

# Override some variables for binmake.gmk.
#
INCFLAGS := -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include

LDLIBS   := $(UILIBS) \
            -lG4intercoms \
            -lG4materials \
            -lG4particles \
            -lG4piiman \
            -lG4track

ifneq (,$(findstring SGI,$(G4SYSTEM)))
        LDLIBS += -lG4particles -lG4global
else
        LDLIBS += -lG4global
endif

include $(G4INSTALL)/config/binmake.gmk
