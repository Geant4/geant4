# $Id: GNUmakefile,v 1.1 1999-01-07 16:13:31 gunter Exp $
# -----------------------------------------------------------
# GNUmakefile for hadronic processes library.  G.Folger 10-Dec-97
# --------------------------------------------------------------

name := G4hadronic_proc

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
