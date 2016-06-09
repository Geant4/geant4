# $Id: GNUmakefile,v 1.2 2006/08/08 05:26:19 kmura Exp $
# $Name: geant4-09-01 $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================
include ../../../../config/config.gmk

LIBRARY  := libezgeom#

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/module.gmk

install:
	@echo ... intall $(lib_name) into $(Q_LIBDIR)
	@if [ ! -d $(Q_LIBDIR) ]; then install -d $(Q_LIBDIR); fi
	@install -m 755 $(LIBRARY).$(soext) $(Q_LIBDIR)

