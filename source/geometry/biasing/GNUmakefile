# $Id: GNUmakefile,v 1.3 2003/10/01 14:51:38 gcosmo Exp $
# --------------------------------------------------------------------
# GNUmakefile for geometry/biasing library.  Gabriele Cosmo, 25/03/02.
# --------------------------------------------------------------------

name := G4geombias

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include \
	    -I$(G4BASE)/materials/include \
	    -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

include $(G4INSTALL)/config/common.gmk
