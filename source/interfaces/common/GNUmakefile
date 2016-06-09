# $Id: GNUmakefile,v 1.5 2010/01/06 14:08:20 lgarnier Exp $
# -------------------------------------------------------------
# GNUmakefile for interfaces/common library.  John Allison, 6/7/98.

name := G4UIcommon

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

# For debug mode
# CPPFLAGS += -DG4DEBUG_INTERFACES_COMMON

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/G4UI_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/interfaces/basic/include

include $(G4INSTALL)/config/common.gmk
