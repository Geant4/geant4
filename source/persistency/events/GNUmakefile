# $Id: GNUmakefile,v 1.1 1999/01/07 16:10:45 gunter Exp $

name := G4pevents
G4ODBMS := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif
include $(G4INSTALL)/config/architecture.gmk

# $(G4TMPDIR) will be overrided in common.gmk later.
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS  += \
             -Iinclude \
             -I$(G4TMPDIR) \
             -I$(G4INSTALL)/source/run/include \
             -I$(G4INSTALL)/source/event/include \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/global/HEPRandom/include \
             -I$(G4INSTALL)/source/global/HEPGeometry/include \
             -I$(G4INSTALL)/source/geometry/management/include \
             -I$(G4INSTALL)/source/geometry/volumes/include \
             -I$(G4INSTALL)/source/geometry/solids/CSG/include \
             -I$(G4INSTALL)/source/events+tracks/include \
             -I$(G4INSTALL)/source/processes/management/include \
             -I$(G4INSTALL)/source/materials/include \
             -I$(G4INSTALL)/source/particles/management/include \
             -I$(G4INSTALL)/source/graphics_reps/include \
             -I$(G4INSTALL)/source/tracking/include \
             -I$(G4INSTALL)/source/track/include \
             -I$(G4INSTALL)/source/persistency/management/include \
             -I$(G4INSTALL)/source/digits+hits/detector/include \
             -I$(G4INSTALL)/source/digits+hits/hits/include \
             -I$(G4INSTALL)/source/digits+hits/digits/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit DDL dependencies:
#  G4PEvent.ddl includes G4PPrimaryVertex.hh, so G4PPrimaryVertex.hh
#  should be created before G4PEvent.hh.

$(G4TMPDIR)/G4PEvent.hh: $(G4TMPDIR)/G4PPrimaryVertex.hh
#--------------------------------------------------------------------

