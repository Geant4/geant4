# $Id: GNUmakefile,v 1.1.1.1 2000/06/21 12:36:24 barrand Exp $
# -------------------------------------------------------------
# GNUmakefile for visualization library.  John Allison, 5/7/95.

name := G4analysis_management

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/analysis.gmk
G4TMPDIR = $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS += -I$(G4BASE)/global/management/include

include $(G4INSTALL)/config/common.gmk

