# $Id: GNUmakefile,v 1.1 2004-03-05 13:29:38 mkossov Exp $
# ---------------------------------------------------------------------
# GNUmakefile for electroweak library. Mikhail Kossov, 22/1/04.
# ---------------------------------------------------------------------

name := G4photoleptohadr_cs

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
