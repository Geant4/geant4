# $Id: GNUmakefile,v 1.2 1999/10/08 15:49:41 morita Exp $

name := G4pcsg
G4ODBMS := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4OODDLX_BOOT := $(G4TMPSCHEMA_BOOT)

# $(G4TMPDIR) will be overriden in common.gmk later.
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS  += \
             -Iinclude \
             -I$(G4TMPDIR) \
             -I$(G4TMP)/$(G4SYSTEM)/G4pgeomn \
             -I../management/include \
             -I$(G4BASE)/persistency/management/include \
             -I$(G4BASE)/run/include \
             -I$(G4BASE)/intercoms/include \
             -I$(G4BASE)/graphics_reps/include \
             -I$(G4BASE)/global/management/include \
             -I$(G4BASE)/global/HEPRandom/include \
             -I$(G4BASE)/global/HEPGeometry/include \
             -I$(G4BASE)/geometry/management/include \
             -I$(G4BASE)/geometry/magneticfield/include \
             -I$(G4BASE)/geometry/volumes/include \
             -I$(G4BASE)/geometry/solids/CSG/include

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

