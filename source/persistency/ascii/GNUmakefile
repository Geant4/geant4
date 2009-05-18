# $Id: GNUmakefile,v 1.5 2009-05-18 15:46:06 arce Exp $
# -----------------------------------------------------------------------
# GNUmakefile for persistency library.  Gabriele Cosmo, 16/11/96.
# -----------------------------------------------------------------------

name := G4geomtext

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/solids/BREPS/include \
            -I$(G4BASE)/geometry/solids/specific/include \
            -I$(G4BASE)/geometry/solids/Boolean/include \
            -I$(G4BASE)/geometry/divisions/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/run/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/scoring/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/navigation/include \
	    -I$(G4BASE)/track/include

include $(G4INSTALL)/config/common.gmk
