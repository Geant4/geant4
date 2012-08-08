# $Id: GNUmakefile,v 1.4 2010-11-08 10:38:44 maire Exp $
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

