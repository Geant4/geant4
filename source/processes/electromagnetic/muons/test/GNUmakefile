# $Id: GNUmakefile,v 1.1 1999-01-08 16:32:22 gunter Exp $
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
INCFLAGS := \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/digits+hits/detector/include \
            -I$(G4BASE)/digits+hits/digits/include \
            -I$(G4BASE)/digits+hits/hits/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/general/include \
            -I$(G4BASE)/processes/electromagnetic/standard/include \
            -I$(G4BASE)/processes/electromagnetic/muons/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/xrays/include

LDLIBS   := \
            -lG4event \
            -lG4tracking \
            -lG4electromagnetic \
            -lG4track \
            -lG4particles \
            -lG4decay \
            -lG4procman \
            -lG4materials \
            -lG4digits+hits \
            -lG4geometry \
            -lG4graphics_reps \
            -lG4intercoms \
            -lG4global

include $(G4INSTALL)/config/binmake.gmk
