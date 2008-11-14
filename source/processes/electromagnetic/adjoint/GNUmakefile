# $Id: GNUmakefile,v 1.1 2008-11-14 19:54:40 gcosmo Exp $
# --------------------------------------------------------------------
# GNUmakefile for electromagnetic sub-library.  G.Cosmo, 14/11/2008.
# --------------------------------------------------------------------

name := G4emadjoint

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/standard/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/adjoint/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk

