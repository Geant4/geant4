# $Id: GNUmakefile,v 1.1 2006/04/25 09:00:01 kmura Exp $
# $Name: geant4-08-01 $
# ===========================================================
#   Makefile for building Geant4Py site-modules
# ===========================================================
include ../../../../config/config.gmk

MODULE  := libezgeom#

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/module.gmk

install:
	@echo ... intall $(lib_name) into $(Q_LIBDIR)
	@if [ ! -d $(Q_LIBDIR) ]; then install -d $(Q_LIBDIR); fi
	@install -m 755 $(MODULE).$(soext) $(Q_LIBDIR)

