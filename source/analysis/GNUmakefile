# $Id: GNUmakefile,v 1.4 2000/12/12 16:29:24 gunter Exp $
# -----------------------------------------------------------------------
# GNUmakefile for visualization and modeling.  John Allison, 5/7/95.
# Modeling is always made.  Others by environment.  John Allison  4/7/98.
# -----------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4analysis

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/analysis.gmk

G4LIBDIR := $(G4LIB)/$(G4SYSTEM)
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(name)

SUBDIRS += management 
SUBLIBS += G4analysis_management

ifdef G4ANALYSIS_BUILD_LIZARD
  SUBDIRS += Lizard
  SUBLIBS += G4Lizard
endif

ifdef G4ANALYSIS_BUILD_JAS
  SUBDIRS += jas
  SUBLIBS += G4jas
endif

ifdef G4ANALYSIS_BUILD_LAB
  SUBDIRS += Lab
  SUBLIBS += G4Lab
endif

.PHONY: granular glob global clean

granular:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done

includes:
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) $@ ); done

glob global:
	@echo "ERROR: The global target not relevant for the analysis category."

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done
