# $Id: GNUmakefile,v 1.8 2001/03/15 22:58:54 morita Exp $

name := G4pgeomBoolean
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
             -I$(G4BASE)/persistency/global/include \
             -I$(G4BASE)/persistency/geometry/global/include \
             -I$(G4BASE)/persistency/management/include \
             -I$(G4BASE)/global/management/include \
             -I$(G4BASE)/global/HEPGeometry/include \
             -I$(G4BASE)/geometry/management/include \
             -I$(G4BASE)/geometry/magneticfield/include \
             -I$(G4BASE)/geometry/volumes/include \
             -I$(G4BASE)/geometry/solids/Boolean/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit DDL dependencies:

#--------------------------------------------------------------------

