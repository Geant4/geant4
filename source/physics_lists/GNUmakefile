# $Id: GNUmakefile,v 1.6 2008-12-05 15:14:20 gcosmo Exp $
# ----------------------------------------------------------------------
# GNUmakefile for physicslists library, Gunter Folger 25-October 2006
# ----------------------------------------------------------------------

name := G4physicslists

SUBDIRS  = lists 
SUBDIRS += builders

GLOBLIBS  = libG4run.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

SUBLIBS  = G4phys_lists 
SUBLIBS += G4phys_builders 

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
