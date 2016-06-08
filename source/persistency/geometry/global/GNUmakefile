# $Id: GNUmakefile,v 1.6 2001/03/15 22:58:54 morita Exp $

name := G4pgeomGlobal
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
             -I$(G4TMPDIR) \
             -I$(G4INSTALL)/source/persistency/global/include \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/geometry/management/include \
             -I$(G4INSTALL)/source/global/HEPGeometry/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit header file dependencies:

#--------------------------------------------------------------------

