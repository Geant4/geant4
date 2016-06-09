# $Id: GNUmakefile,v 1.29 2007/11/02 16:22:05 gcosmo Exp $
# -----------------------------------------------------------------------
# GNUmakefile for persistency library.  Gabriele Cosmo, 16/11/96.
# -----------------------------------------------------------------------

name := G4persistency

SUBDIRS := mctruth
SUBLIBS = G4mctruth

ifdef G4LIB_BUILD_GDML
  SUBDIRS += gdml
  SUBLIBS += G4gdml
endif

GLOBLIBS  = libG4run.lib libG4event.lib
GLOBLIBS += libG4tracking.lib libG4digits_hits.lib libG4track.lib
GLOBLIBS += libG4particles.lib libG4geometry.lib libG4materials.lib
GLOBLIBS += libG4graphics_reps.lib libG4intercoms.lib libG4global.lib


ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
