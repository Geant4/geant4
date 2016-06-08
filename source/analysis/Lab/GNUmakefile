# $Id: GNUmakefile,v 1.1.1.1 2000/11/07 14:31:38 barrand Exp $

name := G4Lab

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/analysis.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include

ifdef G4ANALYSIS_BUILD_OPEN_SCIENTIST
#  CPPFLAGS +=
endif

include $(G4INSTALL)/config/common.gmk
