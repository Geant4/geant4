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
           -I$(G4BASE)/global/management/include \
           -I$(G4BASE)/intercoms/include \
           -I$(G4BASE)/processes/hadronic/models/generator/string_fragmentation/include \
           -I$(G4BASE)/processes/hadronic/models/quark_molecular_dynamics/util/include \
           -I$(G4BASE)/processes/hadronic/models/quark_molecular_dynamics/body/include 

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/management/include/ \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/models/generator/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/processes/hadronic/models/generator/util/include \
            -I$(G4BASE)/processes/hadronic/models/generator/string_common/include \
            -I$(G4BASE)/processes/hadronic/models/generator/string_fragmentation/include \
            -I$(G4BASE)/processes/hadronic/models/quark_molecular_dynamics/util/include \

LDLIBS   := \
         -lG4hadronic_quark_md_body \
         -lG4hadronic_quark_md_util  \
         -lG4globman

include $(G4INSTALL)/config/binmake.gmk
