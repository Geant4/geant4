# $Id: GNUmakefile,v 2.2 1998/10/21 11:13:12 grichine Exp $
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
