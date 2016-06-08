# $Id: GNUmakefile,v 1.6 2001/03/15 22:58:55 morita Exp $

name := G4pglobal
G4ODBMS := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4OODDLX_BOOT := $(G4SCHEMA_BOOT)

# $(G4TMPDIR) will be overrided in common.gmk later.
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS  += \
             -Iinclude \
             -I$(G4TMPDIR) \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/materials/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit DDL dependencies:

#--------------------------------------------------------------------

