# $Id: GNUmakefile,v 1.4 2003-06-02 17:19:25 gcosmo Exp $
# --------------------------------------------------------------
# Makes libraries for each subdomain: 
#         CSG, specific, Boolean, BREPS, STEP, STEPinterface.
# GNUmakefile for solids libraries.  Gabriele Cosmo, 18/9/96.

include $(G4INSTALL)/config/architecture.gmk
SUBDIRS = CSG specific Boolean BREPS
ifdef G4USE_STEP
  SUBDIRS += STEP STEPinterface
endif

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
