# $Id: GNUmakefile,v 1.9 2004/12/08 19:18:45 daquinog Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 08/04/02.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ../../..
endif
 
include $(G4INSTALL)/config/architecture.gmk

SUBDIRS = B01 B02 B03 

.PHONY : all clean clean_libs

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
