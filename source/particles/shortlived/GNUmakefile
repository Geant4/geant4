# $Id: GNUmakefile,v 1.1 1999-01-07 16:10:37 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for shortlived library.  Gabriele Cosmo, 18/9/96.
# --------------------------------------------------------------

name := G4shortlived

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
