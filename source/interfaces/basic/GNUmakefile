# $Id: GNUmakefile,v 1.4 2007/10/29 18:49:45 asaim Exp $
# -------------------------------------------------------------
# GNUmakefile for interfaces/basic library.  John Allison, 6/7/98.

name := G4UIbasic

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

