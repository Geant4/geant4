# $Id: GNUmakefile,v 1.1 2003-10-07 11:27:26 hpw Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables for binmake.gmk.
#
XXINCFLAGS := \
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
            -I$(G4BASE)/processes/hadronic/models/low_energy/include \
            -I$(G4BASE)/processes/hadronic/models/high_energy/include \
            -I$(G4BASE)/processes/hadronic/models/generator/de_excitation/include \
            -I$(G4BASE)/processes/hadronic/models/generator/diffractive_string/include \
            -I$(G4BASE)/processes/hadronic/models/generator/high_energy/include \
            -I$(G4BASE)/processes/hadronic/models/generator/kinetic_model/include \
            -I$(G4BASE)/processes/hadronic/models/generator/management/include \
            -I$(G4BASE)/processes/hadronic/models/generator/pre_equilibrium/include \
            -I$(G4BASE)/processes/hadronic/models/generator/quark_gluon_string/include \
            -I$(G4BASE)/processes/hadronic/models/generator/string_common/include \
            -I$(G4BASE)/processes/hadronic/models/generator/string_fragmentation/include \
            -I$(G4BASE)/processes/hadronic/models/generator/util/include \
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

#-L/afs/cern.ch/rd44/dev/lib/DEC6-AFS/beta03-01

XXLDLIBS   :=  \
	    -lG4event \
	    -lG4materials \
	    -lG4hadronic_diffstring \
	    -lG4hadronic_stringfrag \
	    -lG4hadronic_string_common \
	    -lG4hadronic_util_gen \
	    -lG4hadronic_man_gen \
	    -lG4hadronic_util_gen \
	    -lG4processes \
	    -lG4particles \
	    -lG4materials \
	    -lG4track \
            -lG4tracking \
            -lG4digits+hits \
	    -lG4geometry \
	    -lG4graphics_reps \
	    -lG4intercoms \
	    -lG4global \
	    -lcomplex \
	    -lCLHEP 
	    
LOADLIBS := $(G4BIN)/$(G4SYSTEM)/hbfinit.o -l$(CLHEP_LIB)
LOADLIBS +=    `cernlib packlib` $(FCLIBS)

include $(G4INSTALL)/config/binmake.gmk
	
$(G4TARGET): $(G4BIN)/$(G4SYSTEM)/hbfinit.o


$(G4BIN)/$(G4SYSTEM)/hbfinit.o: hbfinit.f
	$(FC) -c $^ -o $@
