# $Id: GNUmakefile,v 1.3 2007-11-30 14:28:50 lgarnier Exp $
# -------------------------------------------------------------
# GNUmakefile for interfaces/GAG library.  John Allison, 6/7/98.

name := G4UIGAG

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/interfaces/common/include
CPPFLAGS += -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk
