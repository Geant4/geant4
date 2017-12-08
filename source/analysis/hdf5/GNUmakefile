# $Id$
# --------------------------------------------------------------------
# GNUmakefile for analysis/hdf5 library.
# --------------------------------------------------------------------

name := G4hdf5

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

ifeq ($(G4SYSTEM), WIN32-VC)
CPPFLAGS += -I$(G4BASE)/externals/expat/include
endif
ifdef G4LIB_BUILD_EXPAT
CPPFLAGS += -I$(G4BASE)/externals/expat/include
endif

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/analysis/g4tools/include \
            -I$(G4BASE)/analysis/management/include \
            -I$(G4BASE)/analysis/hntools/include

include $(G4INSTALL)/config/common.gmk
