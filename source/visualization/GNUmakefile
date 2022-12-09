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

SUBDIRS := modeling
SUBLIBS := G4modeling

# For Debug mode
# CPPFLAGS += -DG4DEBUG_VIS_OGL
# CPPFLAGS += -DG4DEBUG_VIS_OI
# CPPFLAGS += -DG4DEBUG_VIS_MANAGEMENT

ifdef G4VIS_BUILD
# Drivers that can be used without external libraries...
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
  SUBDIRS += gMocren
  SUBLIBS += G4GMocren
  SUBDIRS += ToolsSG
  SUBLIBS += G4ToolsSG
# Drivers needing external libraries...
ifdef G4VIS_BUILD_OPENGL_DRIVER
  SUBDIRS += OpenGL
  SUBLIBS += G4OpenGL
endif
ifdef G4VIS_BUILD_OI_DRIVER
  SUBDIRS += OpenInventor
  SUBLIBS += G4OpenInventor
endif
ifdef G4VIS_BUILD_QT3D_DRIVER
  SUBDIRS += Qt3D
  SUBLIBS += G4visQt3D
endif
endif #G4VIS_BUILD

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
