# $Id: GNUmakefile,v 1.3 1999/11/24 20:28:06 morita Exp $

name := G4pdigits
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
             -I$(G4INSTALL)/source/persistency/global/include \
             -I$(G4INSTALL)/source/digits+hits/digits/include \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/graphics_reps/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit header file dependencies:

$(G4SCHEMA_HEADER_DIR)/G4PDCofThisEvent.hh: $(G4SCHEMA_HEADER_DIR)/G4PVDigitsCollection.hh
$(G4SCHEMA_HEADER_DIR)/G4PVDigitsCollection.hh: $(G4SCHEMA_HEADER_DIR)/G4PVDigit.hh
#--------------------------------------------------------------------

