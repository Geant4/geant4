# $Id: GNUmakefile,v 1.10 2005-10-04 09:15:47 gcosmo Exp $
# ----------------------------------------------------------------
# GNUmakefile for parameterisations library.
# ----------------------------------------------------------------

name := G4parmodels

SUBDIRS := gflash  
SUBLIBS = G4gflash

GLOBLIBS  = libG4processes.lib libG4digits_hits.lib libG4track.lib
GLOBLIBS += libG4particles.lib libG4geometry.lib libG4materials.lib
GLOBLIBS += libG4graphics_reps.lib libG4intercoms.lib libG4global.lib

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
