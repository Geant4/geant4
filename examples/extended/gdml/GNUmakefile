# $Id: GNUmakefile,v 1.1.1.1 2002-05-31 00:34:42 radoone Exp $
# --------------------------------------------------------------
# Makes libraries and applications: 
#         libgdml sample N03gdml
# GNUmakefile for gdml libraries and examples.  Radovan Chytracek 30/05/2002.

SUBDIRS = libgdml sample N03gdml

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
