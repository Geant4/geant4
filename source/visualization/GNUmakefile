# $Id: GNUmakefile,v 1.16 2004-06-01 15:41:14 gcosmo Exp $
# -----------------------------------------------------------------------
# GNUmakefile for visualization and modeling.  John Allison, 5/7/95.
# Modeling is always made.  Others by environment.  John Allison  4/7/98.
# -----------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4visualization

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS  = libG4run.lib libG4event.lib
GLOBLIBS += libG4tracking.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

SUBDIRS := modeling
SUBLIBS := G4modeling

# Libraries which can be used without external packages or libraries...
ifdef G4VIS_BUILD
  SUBDIRS += management
  SUBLIBS += G4vis_management
  SUBDIRS += FukuiRenderer
  SUBLIBS += G4FR
  SUBDIRS += HepRep
  SUBLIBS += G4visHepRep
  SUBDIRS += RayTracer
  SUBLIBS += G4RayTracer
  SUBDIRS += VRML
  SUBLIBS += G4VRML
  SUBDIRS += Tree
  SUBLIBS += G4Tree
  SUBDIRS += XXX
  SUBLIBS += G4visXXX

ifdef G4VIS_BUILD_OPENGL_DRIVER
  SUBDIRS += OpenGL
  SUBLIBS += G4OpenGL
endif
ifdef G4VIS_BUILD_OI_DRIVER
  SUBDIRS += OpenInventor
  SUBLIBS += G4OpenInventor
endif
endif #G4VIS_BUILD

.PHONY: granular glob global clean

granular:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

includes:
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) $@ ); done

glob global:
	@echo "ERROR: The global target not relevant for the visualization category."

clean clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) $@); done
