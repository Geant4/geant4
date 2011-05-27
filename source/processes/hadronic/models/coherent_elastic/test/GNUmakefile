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

ifeq ($(G4TARGET), G4LEppTest)
   CPPFLAGS+=-I$(G4INSTALL)/tests/tools/include -DMAKEHBOOK -DMAKEBOOK
   EXTRALIBS+=-lG4hadronic_coherent_elastic
   EXTRALIBS+=-L$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM) -lG4TestTool
   EXTRALIBS+=-L/cern/pro/lib -lpacklib
endif

ifeq ($(G4SYSTEM),Linux-g++)
#     EXTRALIBS+=-lf2c -lgcc
#      EXTRALIBS+=-lg2c -lm -lgcc -lc -lgcc
endif

include $(G4INSTALL)/config/binmake.gmk






