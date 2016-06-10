# $Id: GNUmakefile 68844 2013-04-08 07:21:18Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for global.  Gabriele Cosmo, 26/9/96.
# --------------------------------------------------------------

name := G4global

SUBDIRS = management HEPNumerics HEPRandom
SUBLIBS = G4globman G4hepnumerics G4heprandom

GLOBLIBS = 

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/globlib.gmk

includes::
	@for dir in HEPGeometry; do \
	    (cd $$dir && $(MAKE) $@); \
	 done   
