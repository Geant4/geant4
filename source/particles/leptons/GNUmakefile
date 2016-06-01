# $Id: GNUmakefile,v 2.0 1998/07/02 17:25:41 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for leptons library.  Gabriele Cosmo, 18/9/96.
# --------------------------------------------------------------

name := G4leptons

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include

include $(G4INSTALL)/config/common.gmk
