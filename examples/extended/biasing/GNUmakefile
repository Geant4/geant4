# $Id: GNUmakefile,v 1.6 2002-08-29 15:37:24 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 08/04/02.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ../../..
endif
 
include $(G4INSTALL)/config/architecture.gmk

SUBDIRS = B01 B02 B03 B04 B05 B06 B07 B08 B09

.PHONY : all clean clean_libs

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
