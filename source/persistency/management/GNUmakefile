# $Id: GNUmakefile,v 1.18 2001/03/15 22:58:55 morita Exp $

name := G4pmanagement
G4ODBMS := true
G4ODBMS_INSTALLATION := true
PWD := $(shell (dirname `pwd`/dummy | sed -e 's/ //g'))

ifndef G4INSTALL
  G4INSTALL = $(PWD)/../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4OODDLX_BOOT := $(G4SCHEMA_BOOT)

CPPFLAGS  += \
             -I$(PWD)/include \
             -I$(G4SCHEMA_HEADER_DIR) \
             -I$(G4TMP)/$(G4SYSTEM)/G4pdigits \
             -I$(G4TMP)/$(G4SYSTEM)/G4phits \
             -I$(G4TMP)/$(G4SYSTEM)/G4pevents \
             -I$(G4TMP)/$(G4SYSTEM)/G4prun \
             -I$(G4TMP)/$(G4SYSTEM)/G4pgeomGlobal \
             -I$(G4TMP)/$(G4SYSTEM)/G4pcsg \
             -I$(G4TMP)/$(G4SYSTEM)/G4pgeomBoolean \
             -I$(G4TMP)/$(G4SYSTEM)/G4pgeomn \
             -I$(G4INSTALL)/source/persistency/global/include \
             -I$(G4INSTALL)/source/persistency/digits+hits/hits/include \
             -I$(G4INSTALL)/source/persistency/digits+hits/digits/include \
             -I$(G4INSTALL)/source/persistency/events/include \
             -I$(G4INSTALL)/source/persistency/geometry/global/include \
             -I$(G4INSTALL)/source/persistency/geometry/management/include \
             -I$(G4INSTALL)/source/persistency/run/include \
             -I$(G4INSTALL)/source/digits+hits/hits/include \
             -I$(G4INSTALL)/source/digits+hits/digits/include \
             -I$(G4INSTALL)/source/run/include \
             -I$(G4INSTALL)/source/event/include \
             -I$(G4INSTALL)/source/global/HEPGeometry/include \
             -I$(G4INSTALL)/source/geometry/management/include \
             -I$(G4INSTALL)/source/geometry/volumes/include \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/materials/include \
             -I$(G4INSTALL)/source/intercoms/include


#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit DDL dependencies:

#--------------------------------------------------------------------

