# $Id: GNUmakefile,v 1.1 1999-01-07 16:09:50 gunter Exp $
# ---------------------------------------------------------------
# Makes libraries for each subdomain: barions, ions, mesons.
# GNUmakefile for particles libraries.  Gabriele Cosmo, 18/9/96.
# ---------------------------------------------------------------

SUBDIRS = barions ions mesons

all:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
