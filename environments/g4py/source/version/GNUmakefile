# $Id: GNUmakefile,v 1.1 2006-02-27 09:56:05 kmura Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include  ../../config/config.gmk

.PHONY: install

install:
	@echo ... intall headers into $(Q_INCDIR)/version
	@if [ ! -d $(Q_INCDIR)/version ]; then \
	    install -d $(Q_INCDIR)/version; fi
	@install -m 644 *.hh  $(Q_INCDIR)/version



