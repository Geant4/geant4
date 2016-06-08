# $Id: GNUmakefile,v 1.1 1999/01/07 16:10:49 gunter Exp $

name := G4pcsg
G4ODBMS := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif
include $(G4INSTALL)/config/architecture.gmk

# $(G4TMPDIR) will be overriden in common.gmk later.
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS  += \
             -Iinclude \
             -I$(G4TMPDIR) \
             -I$(G4TMP)/$(G4SYSTEM)/G4pgeomn \
             -I../management/include \
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
# Explicit DDL dependencies:

$(G4TMPDIR)/G4PBox.hh:	$(G4TMPDIR)/G4PCSGSolid.hh
$(G4TMPDIR)/G4PCSGSolid.hh:
$(G4TMPDIR)/G4PCons.hh:	$(G4TMPDIR)/G4PCSGSolid.hh
$(G4TMPDIR)/G4PPara.hh:	$(G4TMPDIR)/G4PCSGSolid.hh
$(G4TMPDIR)/G4PSphere.hh:	$(G4TMPDIR)/G4PCSGSolid.hh
$(G4TMPDIR)/G4PTorus.hh:	$(G4TMPDIR)/G4PCSGSolid.hh
$(G4TMPDIR)/G4PTrap.hh:	$(G4TMPDIR)/G4PCSGSolid.hh
$(G4TMPDIR)/G4PTrd.hh:	$(G4TMPDIR)/G4PCSGSolid.hh
$(G4TMPDIR)/G4PTubs.hh:	$(G4TMPDIR)/G4PCSGSolid.hh

#--------------------------------------------------------------------

