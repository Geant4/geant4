# $Id: GNUmakefile,v 1.2 1999-11-22 10:20:50 hpw Exp $
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
            -I$(G4BASE)/processes/hadronic/models/generator/include \
            -I$(G4BASE)/digits+hits/include \
            -I$(G4BASE)/run/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/solids/BREPS/include \
            -I$(G4BASE)/geometry/solids/STEP//include \
	    -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/global/HEPNumerics/include \
	    -I$(G4BASE)/graphics_reps/include \
	    -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/processes/general/include \
            -I$(G4BASE)/processes/management/include \
	    -I$(G4BASE)/processes/parameterisation/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/body/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/nucleus/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/utils/include \
            -I$(G4BASE)/processes/hadronic/models/generator/util/include \
            -I$(G4BASE)/processes/hadronic/models/generator/management/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/stopping/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/electromagnetic/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/digits+hits/detector/include \
            -I$(G4BASE)/digits+hits/hits/include \
            -I$(G4BASE)/digits+hits/digits/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include 

include $(G4INSTALL)/config/binmake.gmk
