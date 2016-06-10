# $Id: GNUmakefile 71741 2013-06-21 17:38:28Z mkelsey $
# -----------------------------------------------------------
# GNUmakefile for hadronic library.  Gabriele Cosmo, 18/9/96.
#
# Map user environment/GMake variables onto preprocessor debugging flags
#
# 20100922 J. Yarba -- Add include directories for pre-compound model
# 20110919 M. Kelsey -- Add flag for final-state clustering
# 20130621 M. Kelsey -- Move CHECK_ECONS flag to runtime parameters
# -----------------------------------------------------------

name := G4hadronic_bert_cascade

ifdef G4CASCADE_COULOMB_DEV
  CPPFLAGS += -DG4CASCADE_COULOMB_DEV=1
endif

ifdef G4CASCADE_DEBUG_CHARGE
  CPPFLAGS += -DG4CASCADE_DEBUG_CHARGE=1
endif

ifdef G4CASCADE_DEBUG_INTERFACE
  CPPFLAGS += -DG4CASCADE_DEBUG_INTERFACE=1
endif

ifdef G4CASCADE_DEBUG_SAMPLER
  CPPFLAGS += -DG4CASCADE_DEBUG_SAMPLER=1
endif

ifdef G4CASCADE_DEBUG_SORT
  CPPFLAGS += -DG4CASCADE_DEBUG_SORT=1
endif

ifdef G4CASCADE_WATCHER_HIGHZ
  CPPFLAGS += -DG4CASCADE_WATCHER_HIGHZ=1
endif

ifdef G4CASCADE_SKIP_ECONS
  CPPFLAGS += -DG4CASCADE_SKIP_ECONS=1
endif

ifdef G4CASCADE_DO_COALESCENCE
  CPPFLAGS += -DG4CASCADE_DO_COALESCENCE=1
endif

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif
include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
	    -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/management/include/ \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
	    -I$(G4BASE)/processes/hadronic/models/management/include \
	    -I$(G4BASE)/processes/hadronic/models/util/include \
            -I$(G4BASE)/processes/hadronic/models/pre_equilibrium/exciton_model/include/ \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/multifragmentation/include/ \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/fermi_breakup/include/ \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/photon_evaporation/include/ \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/management/include/ \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/evaporation/include/ \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/handler/include/ \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/util/include/ \
	    -I$(G4BASE)/processes/hadronic/models/cascade/cascade/include \
	    -I$(G4BASE)/processes/hadronic/models/cascade/utils/include \
	    -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk

