# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
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
