# $Id: GNUmakefile,v 1.2 2006-11-21 07:34:18 gunter Exp $
# ----------------------------------------------------------------------
# GNUmakefile for physicslists library, Gunter Folger 25-October 2006
# ----------------------------------------------------------------------

name := G4physicslists

SUBDIRS  = lists 
SUBDIRS += builders

SUBLIBS  = G4phys_lists 
SUBLIBS += G4phys_builders 

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
