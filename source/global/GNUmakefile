# $Id: GNUmakefile,v 1.2 2000-11-03 16:40:30 gunter Exp $
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
	@for dir in STLInterface HEPGeometry HEPRandom; do \
	    (cd $$dir && $(MAKE) $@); \
	 done   
