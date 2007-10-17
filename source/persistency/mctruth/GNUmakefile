# $Id: GNUmakefile,v 1.1 2007-10-17 09:33:17 gcosmo Exp $
# -----------------------------------------------------------------------
# GNUmakefile for persistency library.  Gabriele Cosmo, 16/11/96.
# -----------------------------------------------------------------------

name := G4mctruth

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4INSTALL)/source/global/management/include \
            -I$(G4INSTALL)/source/global/HEPGeometry/include \
            -I$(G4INSTALL)/source/digits_hits/hits/include \
            -I$(G4INSTALL)/source/digits_hits/digits/include \
            -I$(G4INSTALL)/source/particles/management/include \
            -I$(G4INSTALL)/source/run/include \
            -I$(G4INSTALL)/source/event/include \
            -I$(G4INSTALL)/source/track/include \
            -I$(G4INSTALL)/source/tracking/include \
            -I$(G4INSTALL)/source/intercoms/include \
            -I$(G4INSTALL)/source/geometry/management/include
ifdef G4LIB_USE_HEPMC
  ifndef HEPMC_DIR
	@echo HEPMC_DIR not defined!
  endif
  CPPFLAGS += -I$(HEPMC_DIR)/include
endif

include $(G4INSTALL)/config/common.gmk

