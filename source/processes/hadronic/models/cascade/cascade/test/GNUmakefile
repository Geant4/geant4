# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables for binmake.gmk.
#
INCFLAGS := \
            -I$(G4BASE)/processes/hadronic/models/neutron_hp/include \
            -I$(G4BASE)/digits+hits/include \
            -I$(G4BASE)/events+tracks/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
	    -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/graphics_reps/include \
	    -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/processes/general/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/models/low_energy/include \
            -I$(G4BASE)/processes/hadronic/models/high_energy/include \
            -I$(G4BASE)/processes/hadronic/models/generator/de_excitation/include \
            -I$(G4BASE)/processes/hadronic/models/generator/util/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/stopping/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/electromagnetic/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include 

LDLIBS   := \
	    -lG4hadronic_neu \
	    -lG4event \
	    -lG4global \
	    -lG4materials \
	    -lG4processes \
	    -lG4particles \
	    -lG4hadronic_deex    \
	    -lG4hadronic_man_gen \
	    -lG4hadronic_util_gen \
	    -lG4processes \
	    -lG4track \
            -lG4tracking \
	    -lG4digits+hits \
	    -lG4geometry \
	    -lG4graphics_reps \
	    -lG4intercoms \
	    -lG4global \
            -l$(RWLIB)

include $(G4INSTALL)/config/binmake.gmk

