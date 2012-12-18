# $Id$
# --------------------------------------------------------------
#

name := G4analysis

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/analysis/include
            
include $(G4INSTALL)/config/common.gmk

.PHONY: global

global: lib

