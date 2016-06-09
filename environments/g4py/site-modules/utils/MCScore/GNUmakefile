# $Id: GNUmakefile,v 1.4 2010-12-02 08:34:07 kmura Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include ../../../config/config.gmk

ifeq ($(Q_WITH_PYTHON3), 0)
  install_dir := $(Q_LIBDIR)/g4py

  include $(G4PY_INSTALL)/config/g4py.gmk
  include $(G4PY_INSTALL)/config/sys/$(Q_SYSTEM).gmk
  include $(G4PY_INSTALL)/config/script-install.gmk

else

.PHONY: all install clean 

all:
	@cd python3 && $(MAKE)

install:
	@cd python3 && $(MAKE) install

clean:
	@cd python3 && $(MAKE) clean

endif
