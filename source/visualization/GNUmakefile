# $Id: GNUmakefile,v 1.9 2001-07-27 22:32:49 johna Exp $
# -----------------------------------------------------------------------
# GNUmakefile for visualization and modeling.  John Allison, 5/7/95.
# Modeling is always made.  Others by environment.  John Allison  4/7/98.
# -----------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4visualization

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

G4LIBDIR := $(G4LIB)/$(G4SYSTEM)
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(name)

SUBDIRS := modeling
SUBLIBS := G4modeling

# Libraries which can be used without external packages or libraries...
ifdef G4VIS_BUILD
  SUBDIRS += management
  SUBLIBS += G4vis_management
  SUBDIRS += FukuiRenderer
  SUBLIBS += G4FR
  SUBDIRS += RayTracer
  SUBLIBS += G4RayTracer
  SUBDIRS += VRML
  SUBLIBS += G4VRML
  SUBDIRS += Tree
  SUBLIBS += G4Tree

ifdef G4VIS_BUILD_OPACS_DRIVER
  SUBDIRS += OPACS
  SUBLIBS += G4OPACS
endif
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
