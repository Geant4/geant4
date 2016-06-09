# $Id: GNUmakefile,v 1.3 2008-12-03 07:24:21 kmura Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include  ../../config/config.gmk

.PHONY: install

install:
	@echo ... intall headers into $(Q_INCDIR)/boost
	@if [ ! -d $(Q_INCDIR)/boost ]; then \
	    install -d $(Q_INCDIR)/boost; fi
	@if [ ! -d $(Q_INCDIR)/boost/detail ]; then \
	    install -d $(Q_INCDIR)/boost/detail; fi
	@install -m 644 *.hpp *.hh  $(Q_INCDIR)/boost
	@install -m 644 detail/*.hpp  $(Q_INCDIR)/boost/detail

