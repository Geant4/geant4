# $Id: GNUmakefile,v 1.1 1999-01-07 16:05:17 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ../..
endif
 
include $(G4INSTALL)/config/architecture.gmk

SUBDIRS = A01 A02 A03 A04 A05 A06

.PHONY : all clean clean_libs

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
