# $Id: GNUmakefile 66356 2012-12-18 09:02:32Z gcosmo $
# --------------------------------------------------------------------
# GNUmakefile for management/csv sub-library.
# --------------------------------------------------------------------

name := G4csv

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/analysis/g4tools/include \
            -I$(G4BASE)/analysis/management/include \
            -I$(G4BASE)/analysis/hntools/include

ifdef G4LIB_BUILD_ZLIB
  CPPFLAGS += -I$(G4BASE)/externals/zlib/include
endif

include $(G4INSTALL)/config/common.gmk
