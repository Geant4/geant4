# $Id: GNUmakefile,v 1.2 2000-04-07 13:03:54 gcosmo Exp $
# --------------------------------------------------------------
# Makes libraries for each subdomain: 
#         STEP, BREPS, CSG, STEPinterface, Boolean, specific.
# GNUmakefile for solids libraries.  Gabriele Cosmo, 18/9/96.

SUBDIRS = CSG specific Boolean BREPS STEP STEPinterface

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
