# $Id: GNUmakefile,v 2.5 1998/11/10 18:20:53 morita Exp $

name := G4pgeomn
G4ODBMS := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif
include $(G4INSTALL)/config/architecture.gmk

# $(G4TMPDIR) will be overriden in common.gmk later.
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS  += \
             -Iinclude \
             -I$(G4TMPDIR) \
             -I../volumes/include \
             -I$(G4TMP)/$(G4SYSTEM)/G4pvol \
             -I$(G4INSTALL)/source/persistency/management/include \
             -I$(G4INSTALL)/source/run/include \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/global/HEPRandom/include \
             -I$(G4INSTALL)/source/global/HEPGeometry/include \
             -I$(G4INSTALL)/source/geometry/management/include \
             -I$(G4INSTALL)/source/geometry/magneticfield/include \
             -I$(G4INSTALL)/source/geometry/volumes/include \
             -I$(G4INSTALL)/source/geometry/solids/CSG/include \
             -I$(G4INSTALL)/source/materials/include \
             -I$(G4INSTALL)/source/graphics_reps/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit header file dependencies:

$(G4TMPDIR)/G4PPVParameterised.hh:	$(G4TMPDIR)/G4PPVReplica.hh
$(G4TMPDIR)/G4PPVReplica.hh:	$(G4TMPDIR)/G4PVPhysicalVolume.hh
$(G4TMPDIR)/G4PPVPlacement.hh:	$(G4TMPDIR)/G4PVPhysicalVolume.hh
$(G4TMPDIR)/G4PVPhysicalVolume.hh:	$(G4TMPDIR)/G4PVSolid.hh

#--------------------------------------------------------------------

