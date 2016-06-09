# $Id: GNUmakefile,v 1.3 2006/06/14 08:23:52 kmura Exp $
# $Name: geant4-08-01 $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include ../../../config/config.gmk

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/sys/$(Q_SYSTEM).gmk

install_dir := $(Q_LIBDIR)/EZsim

pys  := $(wildcard *.py)
pycs := $(patsubst %.py, %.pyc, $(pys))
pyos := $(patsubst %.py, %.pyo, $(pys))

.PHONY: clean

all: $(pycs) $(pyos)


install: $(pycs) $(pyos)
	@echo ... intall "*.py" into $(install_dir)
	@if [ ! -d $(install_dir) ]; then install -d $(install_dir); fi
	@install -m 755 $(pys) $(pycs) $(pyos) $(install_dir)

clean:
	@-\rm -f $(pycs) $(pyos)

