# $Id: GNUmakefile,v 1.3 2000/04/07 14:52:49 johna Exp $

name := G4OPACS

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/specific/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/interfaces/common/include
CPPFLAGS += -I$(G4BASE)/intercoms/include

ifdef G4VIS_BUILD_OPACS_DRIVER
  CPPFLAGS += -I$(COROOT)/include -I$(XXROOT)/include -I$(WOROOT)/include\
              -I$(GOROOT)/include -I$(XOROOT)/include
ifeq ($(G4_OPACS_WIDGET_SET),Xm)
  CPPFLAGS += -DHAS_XM  $(X11FLAGS) $(XMFLAGS)  
else
  CPPFLAGS += -DHAS_XAW $(X11FLAGS) $(XAWFLAGS) 
endif
endif

include $(G4INSTALL)/config/common.gmk
