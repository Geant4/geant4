# $Id: GNUmakefile,v 1.1 2006/07/12 02:01:07 kmura Exp $
# $Name: geant4-09-00 $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include ../../../config/config.gmk

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/sys/$(Q_SYSTEM).gmk

install_dir := $(Q_LIBDIR)

pys  := $(wildcard *.py)
pycs := $(patsubst %.py, %.pyc, $(pys))
pyos := $(patsubst %.py, %.pyo, $(pys))

.PHONY: clean

all: $(pycs) $(pyos)


install: $(pycs) $(pyos)
	@echo ... intall "*.py" into $(install_dir)
	@if [ ! -d $(install_dir) ]; then install -d $(install_dir); fi
	@install -m 644 $(pys) $(pycs) $(pyos) $(install_dir)

clean:
	@-\rm -f $(pycs) $(pyos)

