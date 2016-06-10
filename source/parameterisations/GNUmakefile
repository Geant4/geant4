# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
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
