# $Id: GNUmakefile 110061 2018-05-15 09:06:38Z gcosmo $
# -----------------------------------------------------------------------
# GNUmakefile for Geant4 external packages
# Gunter Folger 26-April-2011 
# -----------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4externals

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

SUBDIRS :=
SUBLIBS :=

ifdef G4LIB_USE_CLHEP
  SUBDIRS += clhep
  SUBLIBS += G4clhep
endif
ifdef G4LIB_USE_EXPAT
  SUBDIRS += expat
  SUBLIBS += G4expat
endif
ifdef G4LIB_USE_ZLIB
  SUBDIRS += zlib
  SUBLIBS += G4zlib
endif

.PHONY: granular obj glob global clean

glob global: granular

obj:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) obj); done

granular:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

includes:
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) $@ ); done

clean clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) $@); done
