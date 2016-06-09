# $Id: GNUmakefile,v 1.10 2008/07/02 14:15:16 gcosmo Exp $
# -----------------------------------------------------------------------
# GNUmakefile for persistency library.  Gabriele Cosmo, 16/11/96.
# -----------------------------------------------------------------------

name := G4gdml

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4INSTALL)/source/global/management/include \
            -I$(G4INSTALL)/source/global/HEPGeometry/include \
            -I$(G4INSTALL)/source/global/HEPNumerics/include \
            -I$(G4INSTALL)/source/digits_hits/hits/include \
            -I$(G4INSTALL)/source/digits_hits/digits/include \
            -I$(G4INSTALL)/source/particles/management/include \
            -I$(G4INSTALL)/source/run/include \
            -I$(G4INSTALL)/source/event/include \
            -I$(G4INSTALL)/source/track/include \
            -I$(G4INSTALL)/source/tracking/include \
            -I$(G4INSTALL)/source/intercoms/include \
            -I$(G4INSTALL)/source/geometry/management/include \
	    -I$(G4INSTALL)/source/geometry/volumes/include \
            -I$(G4INSTALL)/source/geometry/navigation/include \
            -I$(G4INSTALL)/source/geometry/solids/CSG/include \
            -I$(G4INSTALL)/source/geometry/solids/Boolean/include \
            -I$(G4INSTALL)/source/geometry/solids/specific/include \
            -I$(G4INSTALL)/source/geometry/divisions/include \
            -I$(G4INSTALL)/source/materials/include \
            -I$(G4INSTALL)/source/graphics_reps/include
 
ifdef G4LIB_USE_GDML
  ifndef XERCESCROOT
	@echo XERCESCROOT not defined!
  endif
endif

include $(G4INSTALL)/config/common.gmk
