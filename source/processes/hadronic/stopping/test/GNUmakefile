# $Id: GNUmakefile,v 1.1 1999-01-08 16:33:46 gunter Exp $
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

ifeq ($(G4TARGET),G4PiMinusAbsorptionAtRestTest)
  EXTRALIBDIR := -L/cern/pro/lib
  EXTRALIBS := -lpacklib $(FCLIBS) 
    ifeq ($(G4SYSTEM),SUN-CC)
      EXTRALIBS := -lpacklib -lw -lgen -lsocket -lnsl -lintl $(FCLIBS) 
    endif
  LDFLAGS += $(EXTRALIBDIR)
  LOADLIBS += $(EXTRALIBS) 
endif

ifeq ($(G4TARGET),G4KaonMinusAbsorptionAtRestTest)
  EXTRALIBDIR := -L/cern/pro/lib
  EXTRALIBS := -lpacklib $(FCLIBS) 
  LDFLAGS += $(EXTRALIBDIR)
  LOADLIBS += $(EXTRALIBS) 
endif

# Override some variables for binmake.gmk.
#
INCFLAGS += \
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
            -I$(G4BASE)/digits+hits/include \
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
            -I$(G4BASE)/processes/hadronic/stopping/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/models/generator/de_excitation/include \
            -I$(G4BASE)/processes/hadronic/models/generator/util/include \
            -I$(G4BASE)/processes/hadronic/management/include

LDLIBS   := \
            -lG4event \
            -lG4tracking \
            -lG4hadronic_stop \
            -lG4hadronic_util \
            -lG4hadronic_deex \
	    -lG4hadronic_util_gen \
            -lG4hadronic_mgt \
            -lG4procman \
            -lG4track \
            -lG4particles \
            -lG4materials \
            -lG4geometry \
            -lG4graphics_reps \
            -lG4intercoms \
            -lG4global

include $(G4INSTALL)/config/binmake.gmk






