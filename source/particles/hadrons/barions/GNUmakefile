# $Id: GNUmakefile,v 1.2 1999-10-07 14:10:23 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for barions library.  Gabriele Cosmo, 18/9/96.
# --------------------------------------------------------------

name := G4baryons

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include

include $(G4INSTALL)/config/common.gmk
