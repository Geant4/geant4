# $Id: GNUmakefile,v 1.1 1999-01-08 16:34:44 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------
# You must have unique G4EXAMPLE_FDID assigned by your local 
# manager of Objectivity/DB lock server.
#
# G4EXAMPLE_BOOT is defined as $HOME/G4EXAMPLE by default.

name := test04
G4ODBMS := true
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all cleandb
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk

  CPPFLAGS += \
              -Iinclude \
              -I$(G4TMPDIR)

test:
	@echo $(CPPFLAGS)

# --------------------------------------------------------------

# override the target all, and warn if G4ODBMS flag is not defined.
ifndef G4ODBMS
all:
	@echo Error - Persistency requires G4ODBMS and related database setup.
	exit 1
endif

cleandb:
ifndef G4EXAMPLE_FDID
	@echo G4EXAMPLE_FDID is not defined.  Stop.
	exit 1
endif
	@rm -f $(G4EXAMPLE_BOOT_DIR)/G4EXAMPLE
	@rm -f $(G4EXAMPLE_BOOT_DIR)/G4EXAMPLE.FDDB
	@rm -f $(G4EXAMPLE_BOOT_DIR)/*.G4EXAMPLE.DB
	@echo Copying schema file to $(G4EXAMPLE_BOOT) for FDID $(G4EXAMPLE_FDID)
	@$(HEP_ODBMS_DIR)/etc/getdb $(G4SCHEMA_DIR)/G4SCHEMA \
	  $(G4EXAMPLE_BOOT) $(G4EXAMPLE_FDID)

