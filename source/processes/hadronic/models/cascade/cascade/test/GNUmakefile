# $Id: GNUmakefile,v 1.3 1999-11-15 12:00:52 miheikki Exp $
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

# Spesific for cascade test
INCFLAGS += \
	-I$(G4BASE)/processes/hadronic/management/include \
	-I$(G4BASE)/processes/hadronic/stopping/include \
	-I$(G4BASE)/processes/hadronic/util/include \
	-I$(G4BASE)/processes/hadronic/processes/include \
	-I$(G4BASE)/processes/hadronic/cross_sections/include \
	-I$(G4BASE)/processes/hadronic/models/neutron_hp/include \
	-I$(G4BASE)/processes/hadronic/models/low_energy/include \
	-I$(G4BASE)/processes/hadronic/models/high_energy/include \
	-I$(G4BASE)/processes/hadronic/models/generator/management/include \
	-I$(G4BASE)/processes/hadronic/models/generator/de_excitation/include \
	-I$(G4BASE)/processes/hadronic/models/generator/diffractive_string/include \
	-I$(G4BASE)/processes/hadronic/models/generator/high_energy/include \
	-I$(G4BASE)/processes/hadronic/models/generator/kinetic_model/include \
	-I$(G4BASE)/processes/hadronic/models/generator/pre_equilibrium/include \
	-I$(G4BASE)/processes/hadronic/models/generator/string_common/include \
	-I$(G4BASE)/processes/hadronic/models/generator/string_fragmentation/include \
	-I$(G4BASE)/processes/hadronic/models/generator/util/include \
	-I$(G4BASE)/processes/hadronic/models/cascade/cascade/include \
	-I$(G4BASE)/processes/hadronic/models/cascade/evaporation/include \

LDLIBS   := \
            -lG4track \
            -lG4partman \
            -lG4intercoms \
            -lG4materials \
            -lG4global

# Spesific for cascade test
LDLIBS   := \
	-lG4hadronic_mgt \
	-lG4hadronic_stop \
	-lG4hadronic_util \
	-lG4hadronic_proc \
	-lG4hadronic_xsect \
	-lG4hadronic_neu \
	-lG4hadronic_LE \
	-lG4hadronic_HE \
	-lG4hadronic_man_gen \
	-lG4hadronic_deex \
	-lG4hadronic_diffstring \
	-lG4hadronic_HE_gen \
	-lG4hadronic_kinetic \
	-lG4hadronic_preequ \
	-lG4hadronic_string_common \
	-lG4hadronic_stringfrag \
	-lG4hadronic_util_gen \
	-lG4hadronic_cascade_bertini \
	-lG4hadronic_cascade_evaporation \

include $(G4INSTALL)/config/binmake.gmk






