# $Id: GNUmakefile,v 1.3 2004-02-08 14:01:50 hpw Exp $
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
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/models/util/include \
            -I$(G4BASE)/processes/hadronic/models/management/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/evaporation/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/fermi_breakup/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/fission/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/gem_evaporation/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/handler/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/management/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/multifragmentation/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/photon_evaporation/include \
              -I$(G4BASE)/processes/hadronic/models/de_excitation/util/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/models/generator/de_excitation/include \
            -I$(G4BASE)/processes/hadronic/models/generator/util/include \
            -I$(G4BASE)/processes/hadronic/management/include


include $(G4INSTALL)/config/binmake.gmk






