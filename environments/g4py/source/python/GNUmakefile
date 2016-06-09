# $Id: GNUmakefile,v 1.3 2006/06/29 01:39:00 kmura Exp $
# $Name: geant4-08-02 $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
include ../../config/config.gmk

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/sys/$(Q_SYSTEM).gmk

G4PY_LIBDIR  := $(Q_LIBDIR)/Geant4

pys  := $(wildcard *.py)
pycs := $(patsubst %.py, %.pyc, $(pys))
pyos := $(patsubst %.py, %.pyo, $(pys))

.PHONY: clean

all: $(pycs) $(pyos)

install: $(pycs) $(pyos)
	@echo ... intall "*.py" into $(G4PY_LIBDIR)
	@if [ ! -d $(G4PY_LIBDIR) ]; then install -d $(G4PY_LIBDIR); fi
	@install -m 644 $(pys) $(pycs) $(pyos) $(G4PY_LIBDIR)

clean:
	@-\rm -f $(pycs) $(pyos)

