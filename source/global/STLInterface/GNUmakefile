# $Id: GNUmakefile,v 1.1 2000/11/03 16:39:29 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for STLInterface.  Gunter Folger 3-Nov-2000.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4DUMMY_VARIABLE:=$(shell [ ! -d $(G4INCLUDE)/g4rw ] && mkdir -p $(G4INCLUDE)/g4rw )

installed_includes:=$(foreach file,$(wildcard g4rw/*),$(shell test -f $(file) && echo $(file)))

installed_includes:=$(foreach file,$(installed_includes),$(G4INCLUDE)/$(file))

includes:	$(installed_includes)

$(installed_includes): $(G4INCLUDE)/% : %
	@cp -p $< $@
	
