# $Id: GNUmakefile,v 1.9 1999/12/01 15:21:05 morita Exp $

name := G4pevents
G4ODBMS := true
G4ODBMS_INSTALLATION := true


ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4OODDLX_BOOT := $(G4TMPSCHEMA_BOOT)

CPPFLAGS  += \
             -Iinclude \
             -I$(G4SCHEMA_HEADER_DIR) \
             -I$(G4TMP)/$(G4SYSTEM)/G4phits \
             -I$(G4TMP)/$(G4SYSTEM)/G4pdigits \
             -I$(G4INSTALL)/source/persistency/global/include \
             -I$(G4INSTALL)/source/persistency/digits+hits/digits/include \
             -I$(G4INSTALL)/source/persistency/digits+hits/hits/include \
             -I$(G4INSTALL)/source/persistency/management/include \
             -I$(G4INSTALL)/source/event/include \
             -I$(G4INSTALL)/source/tracking/include \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/graphics_reps/include \
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

$(G4SCHEMA_HEADER_DIR)/G4PEvent.hh: $(G4SCHEMA_HEADER_DIR)/G4PPrimaryVertex.hh
#--------------------------------------------------------------------

