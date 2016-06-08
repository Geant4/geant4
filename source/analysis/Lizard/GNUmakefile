# $Id: GNUmakefile,v 1.1.1.1 2000/10/31 21:42:25 barrand Exp $

name := G4Lizard

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/analysis.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include

ifdef G4ANALYSIS_BUILD_LIZARD
#  CPPFLAGS +=
endif

include $(G4INSTALL)/config/common.gmk
