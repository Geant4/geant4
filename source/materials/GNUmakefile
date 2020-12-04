# --------------------------------------------------------------
# GNUmakefile for materials library.  Gabriele Cosmo, 18/9/96.
# --------------------------------------------------------------

name := G4materials


ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4intercoms.lib libG4global.lib

ifdef G4LIB_BUILD_ZLIB
  GLOBLIBS += libG4zlib.lib
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/intercoms/include

ifdef G4LIB_BUILD_ZLIB
  CPPFLAGS += -I$(G4BASE)/externals/zlib/include
endif


include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

