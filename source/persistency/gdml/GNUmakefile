# $Id: GNUmakefile 85532 2014-10-30 16:19:11Z gcosmo $
# -----------------------------------------------------------------------
# GNUmakefile for persistency library.  Gabriele Cosmo, 16/11/96.
# -----------------------------------------------------------------------

name := G4gdml

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/global/HEPNumerics/include \
	    -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/digits_hits/hits/include \
            -I$(G4BASE)/digits_hits/digits/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/run/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/intercoms/include \
	    -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/geometry/management/include \
	    -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/solids/Boolean/include \
            -I$(G4BASE)/geometry/solids/specific/include \
            -I$(G4BASE)/geometry/divisions/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/graphics_reps/include
 
ifdef G4LIB_USE_GDML
  ifndef XERCESCROOT
	@echo XERCESCROOT not defined!
  endif
endif

include $(G4INSTALL)/config/common.gmk
