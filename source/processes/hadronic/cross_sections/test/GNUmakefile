# $Id: GNUmakefile,v 1.2 1999-10-06 13:40:38 hpw Exp $
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

ifeq ($(G4TARGET),G4HadronAtRestTest)
  LDSTUBS :=  $(G4TMP)/$(G4SYSTEM)/$(G4TARGET)/G4HadronAtRestTestStub.o
  EXTRALIBDIR := -L/cern/pro/lib
  EXTRALIBS := -lpacklib $(FCLIBS) $(LDSTUBS)
  LFLAGS += $(EXTRALIBDIR)
  LOADLIBS += $(EXTRALIBS)
endif

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
            -I$(G4BASE)/particle+matter/general/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/processes/general/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/electromagnetic/standard/include \
            -I$(G4BASE)/processes/electromagnetic/muons/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/xrays/include \
            -I$(G4BASE)/run/include \
            -I$(G4BASE)/digits+hits/include \
            -I$(G4BASE)/intercoms/include

LDLIBS   := \
	    -lG4hadronic_xsect \
            -lG4processes \
            -lG4track \
            -lG4tracking \
            -lG4particles \
            -lG4materials \
            -lG4geometry \
            -lG4graphics_reps \
            -lG4intercoms \
            -lG4event \
            -lG4digits+hits \
            -lG4global \

include $(G4INSTALL)/config/binmake.gmk
