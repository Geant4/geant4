# $Id: GNUmakefile,v 1.1 1999-01-07 16:14:07 gunter Exp $
# ---------------------------------------------------------------------
# GNUmakefile for photolepton_hadron library.  Gabriele Cosmo, 18/9/96.
# ---------------------------------------------------------------------

name := G4photolepton_hadron

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -g
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
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/models/high_energy/include \
            -I$(G4BASE)/processes/hadronic/models/low_energy/include

include $(G4INSTALL)/config/common.gmk
