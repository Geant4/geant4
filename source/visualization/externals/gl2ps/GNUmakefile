# $Id: GNUmakefile,v 1.5 2009-11-18 17:15:11 gcosmo Exp $
# -------------------------------------------------------------
# GNUmakefile for gl2ps.  Laurent Garnier, 6/2/09.

name := G4gl2ps

# For debug mode
# CPPFLAGS += -DG4DEBUG_VIS_GL2PS

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

GLOBLIBS = libG4zlib.lib

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

ifdef G4LIB_BUILD_ZLIB
  CPPFLAGS += -I$(G4BASE)/visualization/externals/zlib/include
endif

include $(G4INSTALL)/config/common.gmk
