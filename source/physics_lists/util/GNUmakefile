# $Id: GNUmakefile 66704 2013-01-10 18:20:17Z gunter $
# ---------------------------------------------------------------------------
# GNUmakefile for physics_lists/lists library.  Gunter Folger 25-Oct-2006.
# ---------------------------------------------------------------------------

name := G4physlist_util

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4PROCESSES_EXPORT

ifeq ($(G4SYSTEM), WIN32-VC)
CPPFLAGS += -I$(G4BASE)/externals/expat/include
endif

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/models/management/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
