# $Id: GNUmakefile,v 1.7 2000-11-03 16:41:10 gunter Exp $
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

ifdef G4VIS_BUILD
  SUBDIRS += management
  SUBLIBS += G4vis_management
endif

# For DAWN and DAWNFILE drivers
ifdef G4VIS_BUILD_DAWN_DRIVER
  SUBDIRS += FukuiRenderer
  SUBLIBS += G4FR
endif
ifdef G4VIS_BUILD_DAWNFILE_DRIVER
  ifndef   G4VIS_BUILD_DAWN_DRIVER
    SUBDIRS += FukuiRenderer
    SUBLIBS += G4FR
  endif
endif
# End of DAWN and DAWNFILE drivers 

ifdef G4VIS_BUILD_RAYTRACER_DRIVER
  SUBDIRS += RayTracer
  SUBLIBS += G4RayTracer
endif

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

# For VRML and VRMLFILE drivers
ifdef      G4VIS_BUILD_VRML_DRIVER
  SUBDIRS += VRML
  SUBLIBS += G4VRML
endif
ifdef      G4VIS_BUILD_VRMLFILE_DRIVER
  ifndef   G4VIS_BUILD_VRML_DRIVER
    SUBDIRS += VRML
    SUBLIBS += G4VRML
  endif
endif
# End of VRML and VRMLFILE drivers

.PHONY: granular glob global clean

granular:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

includes:
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) $@ ); done

glob global:
	@echo "ERROR: The global target not relevant for the visualization category."

clean clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) $@); done
	
