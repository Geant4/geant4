# $Id: GNUmakefile,v 1.3 2003/10/03 14:34:52 gcosmo Exp $
# ---------------------------------------------------------------
# GNUmakefile for parameterisations library.  G.Cosmo, 14/11/2000
# ---------------------------------------------------------------

name := G4parmodels

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/processes/parameterisation/include \
            -I$(G4BASE)/processes/electromagnetic/standard/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/digits_hits/detector/include \
            -I$(G4BASE)/digits_hits/hits/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global

global: lib
