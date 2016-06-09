# $Id: GNUmakefile,v 1.2 2005/01/26 08:50:58 gcosmo Exp $
# -------------------------------------------------------------
# GNUmakefile for interfaces/GAG library.  John Allison, 6/7/98.

name := G4UIGAG

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/interfaces/common/include
CPPFLAGS += -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk
