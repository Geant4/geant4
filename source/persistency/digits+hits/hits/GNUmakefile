# $Id: GNUmakefile,v 1.1 1999/01/07 16:10:44 gunter Exp $

name := G4phits
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
             -I$(G4INSTALL)/source/persistency/management/include \
             -I$(G4INSTALL)/source/run/include \
             -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/digits+hits/hits/include \
             -I$(G4INSTALL)/source/graphics_reps/include

#--------------------------------------------------------------------

include $(G4INSTALL)/config/G4ODBMS_BUILD.gmk

include $(G4INSTALL)/config/common.gmk

#--------------------------------------------------------------------
# Explicit header file dependencies:

#--------------------------------------------------------------------

