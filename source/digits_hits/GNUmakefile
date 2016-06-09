# $Id: GNUmakefile,v 1.2 2003/10/03 10:24:01 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for digits+hits library.  Makoto Asai, 1/11/96.
# --------------------------------------------------------------

name := G4digits_hits

SUBDIRS := detector hits digits
SUBLIBS = G4detector G4digits G4hits

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk
