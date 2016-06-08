# $Id: GNUmakefile,v 1.10 1999/12/01 15:21:06 morita Exp $

name := G4pgeomn
G4ODBMS := true
G4ODBMS_INSTALLATION := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4OODDLX_BOOT := $(G4TMPSCHEMA_BOOT)

CPPFLAGS  += \
             -Iinclude \
             -I$(G4SCHEMA_HEADER_DIR) \
             -I$(G4TMP)/$(G4SYSTEM)/G4pgeomGlobal \
             -I$(G4TMP)/$(G4SYSTEM)/G4pcsg \
             -I$(G4TMP)/$(G4SYSTEM)/G4pgeomBoolean \
             -I$(G4INSTALL)/source/persistency/global/include \
             -I$(G4INSTALL)/source/persistency/geometry/global/include \
             -I$(G4INSTALL)/source/persistency/management/include \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/global/HEPRandom/include \
             -I$(G4INSTALL)/source/global/HEPGeometry/include \
             -I$(G4INSTALL)/source/geometry/management/include \
             -I$(G4INSTALL)/source/geometry/magneticfield/include \
             -I$(G4INSTALL)/source/geometry/volumes/include \
             -I$(G4INSTALL)/source/geometry/solids/CSG/include \
             -I$(G4INSTALL)/source/geometry/solids/Boolean/include \
             -I$(G4INSTALL)/source/materials/include \
             -I$(G4INSTALL)/source/graphics_reps/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit header file dependencies:

$(G4SCHEMA_HEADER_DIR)/G4PGeometryObjectMap.hh:	\
	$(G4SCHEMA_HEADER_DIR)/G4PVPhysicalVolume.hh \
	$(G4SCHEMA_HEADER_DIR)/G4PLogicalVolume.hh
$(G4SCHEMA_HEADER_DIR)/G4PPVParameterised.hh: \
	$(G4SCHEMA_HEADER_DIR)/G4PPVReplica.hh
$(G4SCHEMA_HEADER_DIR)/G4PPVReplica.hh:	\
	$(G4SCHEMA_HEADER_DIR)/G4PVPhysicalVolume.hh
$(G4SCHEMA_HEADER_DIR)/G4PPVPlacement.hh: \
	$(G4SCHEMA_HEADER_DIR)/G4PVPhysicalVolume.hh
$(G4SCHEMA_HEADER_DIR)/G4PLogicalVolume.hh: \
	$(G4SCHEMA_HEADER_DIR)/G4PVPhysicalVolume.hh
#--------------------------------------------------------------------

