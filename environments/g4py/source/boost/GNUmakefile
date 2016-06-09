# $Id: GNUmakefile,v 1.2 2006/04/25 08:09:43 kmura Exp $
# $Name: geant4-09-00 $
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

