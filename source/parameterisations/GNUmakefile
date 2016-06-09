# $Id: GNUmakefile,v 1.9 2004/12/15 08:49:47 gcosmo Exp $
# ----------------------------------------------------------------
# GNUmakefile for parameterisations library.  HPW 09-Nov-04.
# ----------------------------------------------------------------

name := G4parmodels

SUBDIRS := gflash trd_clusters  
SUBLIBS = G4gflash G4trdmodels

GLOBLIBS  = libG4processes.lib libG4digits_hits.lib libG4track.lib
GLOBLIBS += libG4particles.lib libG4geometry.lib libG4materials.lib
GLOBLIBS += libG4graphics_reps.lib libG4intercoms.lib libG4global.lib

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
