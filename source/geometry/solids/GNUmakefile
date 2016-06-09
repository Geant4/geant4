# $Id: GNUmakefile,v 1.6 2003-10-01 14:51:39 gcosmo Exp $
# --------------------------------------------------------------
# Makes libraries for each subdomain: 
#         CSG, specific, Boolean, BREPS, STEP, STEPinterface.
# GNUmakefile for solids libraries.  Gabriele Cosmo, 18/9/96.

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
SUBDIRS = CSG specific Boolean BREPS

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
