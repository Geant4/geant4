# $Id: GNUmakefile 85530 2014-10-30 16:17:43Z gcosmo $
# --------------------------------------------------------------
# Makes libraries for each subdomain: 
#         CSG, specific, Boolean.
# GNUmakefile for solids libraries.  Gabriele Cosmo, 18/9/96.

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
SUBDIRS = CSG specific Boolean

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
