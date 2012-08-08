# $Id: GNUmakefile,v 1.8 2009/06/18 16:18:41 mantero Exp $
# --------------------------------------------------------------------
# GNUmakefile for electromagnetic sub-library.  John Allison, 25/6/98.
# --------------------------------------------------------------------

name := G4emdna-moltypes

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4EM_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/processes/electromagnetic/dna/management/include \
            -I$(G4BASE)/processes/electromagnetic/dna/molecules/management/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk


