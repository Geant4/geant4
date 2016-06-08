# $Id: GNUmakefile,v 1.2 2001/03/15 22:58:54 morita Exp $

name := G4pspecsolids
G4ODBMS := true
G4ODBMS_INSTALLATION := true
PWD := $(shell (dirname `pwd`/dummy | sed -e 's/ //g'))

ifndef G4INSTALL
  G4INSTALL = $(PWD)/../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4OODDLX_BOOT := $(G4SCHEMA_BOOT)

CPPFLAGS  += \
             -I$(PWD)/include \
             -I$(G4SCHEMA_HEADER_DIR) \
             -I$(G4TMP)/$(G4SYSTEM)/G4pgeomGlobal \
             -I$(G4TMP)/$(G4SYSTEM)/G4pcsg \
             -I$(G4BASE)/persistency/global/include \
             -I$(G4BASE)/persistency/geometry/global/include \
             -I$(G4BASE)/graphics_reps/include \
             -I$(G4BASE)/global/management/include \
             -I$(G4BASE)/global/HEPRandom/include \
             -I$(G4BASE)/global/HEPGeometry/include \
             -I$(G4BASE)/geometry/management/include \
             -I$(G4BASE)/geometry/magneticfield/include \
             -I$(G4BASE)/geometry/volumes/include \
             -I$(G4BASE)/geometry/solids/CSG/include \
             -I$(G4BASE)/geometry/solids/specific/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit DDL dependencies:

##$(G4SCHEMA_HEADER_DIR)/G4PBox.hh:	$(G4SCHEMA_HEADER_DIR)/G4PCSGSolid.hh

#--------------------------------------------------------------------

