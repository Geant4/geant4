# $Id: GNUmakefile,v 1.2 2006/04/25 08:58:18 kmura Exp $
# $Name: geant4-08-01 $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include ../../../config/config.gmk

# -----------------------------------------------------------
# library
# -----------------------------------------------------------
SUBDIR := ezgeom

.PHONY: all cleanall

all:
	@cd $(SUBDIR) && $(MAKE)
	@$(MAKE) $(module_name)

installall:
	@cd $(SUBDIR) && $(MAKE) install
	@$(MAKE) install

cleanall:
	@cd $(SUBDIR) && $(MAKE) clean
	@-rm -f $(module_name)
	@-rm -f *.o

# -----------------------------------------------------------
# module
# -----------------------------------------------------------
# python module name
PACKAGE := EZsim#
MODULE  := EZgeom#

# location of exposured library
G4PY_EXLIB_LIBDIR  := $(Q_LIBDIR)

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/module.gmk
include $(G4PY_INSTALL)/config/site-install.gmk

CXXFLAGS += -I./$(SUBDIR)
LOPT     += -L./$(SUBDIR) -lezgeom

