# $Id: GNUmakefile,v 1.21 2002-11-24 13:45:22 morita Exp $
# -----------------------------------------------------------------------
# GNUmakefile for persistency library.  Gabriele Cosmo, 16/11/96.
# -----------------------------------------------------------------------

name := G4persistency

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4INSTALL)/source/global/management/include \
            -I$(G4INSTALL)/source/global/HEPGeometry/include \
            -I$(G4INSTALL)/source/digits+hits/hits/include \
            -I$(G4INSTALL)/source/digits+hits/digits/include \
            -I$(G4INSTALL)/source/particles/management/include \
            -I$(G4INSTALL)/source/run/include \
            -I$(G4INSTALL)/source/event/include \
            -I$(G4INSTALL)/source/track/include \
            -I$(G4INSTALL)/source/tracking/include \
            -I$(G4INSTALL)/source/intercoms/include \
            -I$(G4INSTALL)/source/geometry/management/include

include $(G4INSTALL)/config/common.gmk
