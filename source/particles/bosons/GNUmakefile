# $Id: GNUmakefile,v 1.4 2010-10-27 07:47:04 gcosmo Exp $
# ----------------------------------------------------------
# GNUmakefile for bosons library.  Gabriele Cosmo, 18/9/96.
# ----------------------------------------------------------

name := G4bosons

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4PARTICLES_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include

include $(G4INSTALL)/config/common.gmk
