# $Id: GNUmakefile,v 1.1 1999-01-07 16:07:22 gunter Exp $
# --------------------------------------------------------------
# Makes libraries for each subdomain: 
#         STEP, BREPS, CSG, STEPinterface.
# GNUmakefile for solids libraries.  Gabriele Cosmo, 18/9/96.

SUBDIRS = CSG Boolean BREPS STEP STEPinterface

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
