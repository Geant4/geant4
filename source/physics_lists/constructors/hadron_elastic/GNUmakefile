# $Id: GNUmakefile 104019 2017-05-08 07:34:08Z gcosmo $
# ---------------------------------------------------------------------------
# GNUmakefile for physics_lists/constructors/electromagnetic library.  
#  Gunter Folger 10-Jan-2012.
# ---------------------------------------------------------------------------

name := G4phys_ctor_helastic

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4PROCESSES_EXPORT

ifeq ($(G4SYSTEM), WIN32-VC)
CPPFLAGS += -I$(G4BASE)/externals/expat/include
endif
ifdef G4LIB_BUILD_EXPAT
CPPFLAGS += -I$(G4BASE)/externals/expat/include
endif

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/run/include \
	    -I$(G4BASE)/physics_lists/builders/include \
	    -I$(G4BASE)/physics_lists/util/include \
	    -I$(G4BASE)/physics_lists/constructors/factory/include \
	    -I$(G4BASE)/physics_lists/constructors/hadron_inelastic/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/transportation/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/stopping/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/models/management/include \
            -I$(G4BASE)/processes/hadronic/models/util/include \
            -I$(G4BASE)/processes/hadronic/models/coherent_elastic/include \
            -I$(G4BASE)/processes/hadronic/models/theo_high_energy/include \
            -I$(G4BASE)/processes/hadronic/models/lend/include \
            -I$(G4BASE)/processes/hadronic/models/particle_hp/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
