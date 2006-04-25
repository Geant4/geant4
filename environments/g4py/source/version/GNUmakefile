# $Id: GNUmakefile,v 1.2 2006-04-25 08:09:45 kmura Exp $
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



