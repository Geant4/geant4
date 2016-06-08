# $Id: GNUmakefile,v 1.11 2001/03/15 22:58:54 morita Exp $

name := G4phits
G4ODBMS := true
G4ODBMS_INSTALLATION := true
PWD := $(shell (dirname `pwd`/dummy | sed -e 's/ //g'))

ifndef G4INSTALL
  G4INSTALL = $(PWD)/../../../..
endif
include $(G4INSTALL)/config/architecture.gmk

G4OODDLX_BOOT := $(G4SCHEMA_BOOT)

CPPFLAGS  += \
             -I$(PWD)/include \
             -I$(G4SCHEMA_HEADER_DIR) \
             -I$(G4INSTALL)/source/persistency/global/include \
             -I$(G4INSTALL)/source/digits+hits/hits/include \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/graphics_reps/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit header file dependencies:

$(G4SCHEMA_HEADER_DIR)/G4PHCofThisEvent.hh: $(G4SCHEMA_HEADER_DIR)/G4PVHitsCollection.hh
$(G4SCHEMA_HEADER_DIR)/G4PVHitsCollection.hh: $(G4SCHEMA_HEADER_DIR)/G4PVHit.hh
#--------------------------------------------------------------------

