# $Id: GNUmakefile,v 1.1.1.1 2002-05-31 00:34:42 radoone Exp $
# --------------------------------------------------------------
# Makes libraries for each GDML subdomain: 
#         saxana schema management processes subscribers.
# GNUmakefile for GDML libraries.  Radovan Chytracek, 30/05/2002.

SUBDIRS = saxana schema management processes subscribers

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
