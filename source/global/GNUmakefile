# $Id: GNUmakefile,v 1.3 2001-10-11 14:04:04 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for global.  Gabriele Cosmo, 26/9/96.
# --------------------------------------------------------------

name := G4global

SUBDIRS = management HEPNumerics
SUBLIBS = G4globman G4hepnumerics

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk

includes::
	@for dir in HEPGeometry HEPRandom; do \
	    (cd $$dir && $(MAKE) $@); \
	 done   
