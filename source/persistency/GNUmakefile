# --------------------------------------------------------------
# GNUmakefile for persistency.  Youhei Morita, 16/06/98.
# --------------------------------------------------------------
# object files are compiled from DDL files as well as from CC files
#
# include/*.ddl -> $(G4TMP)/$(G4SYSTEM)/*.hh     -> $(G4SCHEMA_INCLUDE)/*.hh
#                  $(G4TMP)/$(G4SYSTEM)/*_ref.hh -> $(G4SCHEMA_INCLUDE)/*.hh
#                  $(G4TMP)/$(G4SYSTEM)/*_ddl.cc -> $(G4TMP)/$(G4SYSTEM)/*_ddl.o
# (schema) ------> $(G4TMPSCHEMA_BOOT) -----------> $(G4SCHEMA_BOOT)
#  src/*.cc --------------------------------------> $(G4TMP)/$(G4SYSTEM)/*.o
#
# You must have unique G4SCHEMA_FDID assigned by your local 
# manager of Objectivity/DB lock server.
#

name := G4persistency
G4ODBMS := true
G4ODBMS_INSTALLATION := true

SUBDIRS = global \
          digits+hits/digits digits+hits/hits events run geometry/global \
          geometry/solids/CSG geometry/solids/Boolean geometry/management \
          management
SUBLIBS = G4pglobal G4pdigits G4phits G4pevents G4prun G4pgeomGlobal \
          G4pcsg G4pgeomBoolean G4pgeomn \
          G4pmanagement

ifndef G4INSTALL
  G4INSTALL = ../..
endif
include $(G4INSTALL)/config/globlib.gmk

# --------------------------------------------------------------

.PHONY : all global clean clean_libs cleandb

# if G4ODBMS flag is not defined, override the target "all" then exit 
ifndef G4ODBMS
all:
	@echo Error - Persistency requires G4ODBMS and related database setup.
	exit 1
endif

cleandb:
ifndef G4SCHEMA_FDID
	@echo G4SCHEMA_FDID is not defined.  Stop.
	exit 1
endif
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done
	@rm -f $(G4TMPSCHEMA_BOOT_DIR)/G4TMPSCHEMA
	@rm -f $(G4TMPSCHEMA_BOOT_DIR)/G4TMPSCHEMA.FDDB
	@rm -f $(G4TMPSCHEMA_BOOT_DIR)/*.G4TMPSCHEMA.DB
	@echo Copying schema to $(G4TMPSCHEMA_BOOT) for FDID $(G4SCHEMA_FDID)
	@$(HEP_ODBMS_DIR)/etc/getdb $(HEP_ODBMS_DIR)/schema/HEP_BASE \
	  $(G4TMPSCHEMA_BOOT) $(G4SCHEMA_FDID)

