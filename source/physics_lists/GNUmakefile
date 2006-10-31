# $Id: GNUmakefile,v 1.1 2006-10-31 11:34:03 gunter Exp $
# ----------------------------------------------------------------------
# GNUmakefile for physicslists library, Gjunter Folger 25-October 2006
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
