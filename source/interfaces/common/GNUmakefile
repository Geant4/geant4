# $Id: GNUmakefile,v 1.1 1999-01-07 16:09:36 gunter Exp $
# -------------------------------------------------------------
# GNUmakefile for interfaces/common library.  John Allison, 6/7/98.

name := G4UIcommon

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4UI_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

G4TMPDIR = $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk
