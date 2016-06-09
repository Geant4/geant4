# $Id: GNUmakefile,v 1.3 2006/08/10 08:50:18 kmura Exp $
# $Name: geant4-09-01 $
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

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/module.gmk
include $(G4PY_INSTALL)/config/site-install.gmk

CXXFLAGS += -I./$(SUBDIR)
LOPT     += -L./$(SUBDIR) -lezgeom

