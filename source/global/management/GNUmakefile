# $Id: GNUmakefile,v 1.2 2000-11-03 16:40:30 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for global management.  Gabriele Cosmo, 26/9/96.
# --------------------------------------------------------------

name := G4globman

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/common.gmk


G4DUMMY_VARIABLE:=$(shell [ ! -d $(G4INCLUDE)/g4std ] && mkdir -p $(G4INCLUDE)/g4std )

installed_includes:=$(foreach file,$(wildcard include/g4std/*),$(shell test -f $(file) && echo $(file)))

installed_includes:=$(patsubst include/%,$(G4INCLUDE)/%,$(installed_includes))

includes::	$(installed_includes)

$(installed_includes): $(G4INCLUDE)/% : include/%
	@cp -p $< $@
