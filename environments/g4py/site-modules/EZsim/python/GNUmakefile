# $Id: GNUmakefile,v 1.2 2006-06-13 16:31:17 kmura Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include ../../../config/config.gmk

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/sys/$(Q_SYSTEM).gmk

#install_dir := $(Q_LIBDIR)/site-modules/EZsim
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

