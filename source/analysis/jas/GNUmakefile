# $Id: GNUmakefile,v 1.2 2000/09/14 12:42:32 barrand Exp $

name := G4jas

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/analysis.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include

ifdef G4ANALYSIS_BUILD_JAS
#  CPPFLAGS +=
endif

include $(G4INSTALL)/config/common.gmk
