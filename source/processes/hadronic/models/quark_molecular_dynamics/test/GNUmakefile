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
           -I$(G4BASE)/processes/hadronic/models/quark_molecular_dynamics/util/include \
           -I$(G4BASE)/processes/hadronic/models/quark_molecular_dynamics/body/include 
                      

LDLIBS   := \
         -lG4hadronic_quark_md_body \
         -lG4hadronic_quark_md_util  \
         -lG4globman

include $(G4INSTALL)/config/binmake.gmk
