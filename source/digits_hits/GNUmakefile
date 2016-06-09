# $Id: GNUmakefile,v 1.3 2004/06/01 15:34:51 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for digits+hits library.  Makoto Asai, 1/11/96.
# --------------------------------------------------------------

name := G4digits_hits

SUBDIRS := detector hits digits
SUBLIBS = G4detector G4digits G4hits

GLOBLIBS  = libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4intercoms.lib libG4global.lib

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
