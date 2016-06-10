# $Id: GNUmakefile 66356 2012-12-18 09:02:32Z gcosmo $
# --------------------------------------------------------------------
# GNUmakefile for analysis/management sub-library.
# --------------------------------------------------------------------

name := G4analysismng

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/analysis/g4tools/include

include $(G4INSTALL)/config/common.gmk
