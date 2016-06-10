# $Id:$
# ------------------------------------------------------------
# GNUmakefile for internal tpmalloc library.
# ------------------------------------------------------------

name := G4tpmalloc

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

GLOBLIBS = libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

