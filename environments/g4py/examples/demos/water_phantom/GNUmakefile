# $Id: GNUmakefile,v 1.1 2006-02-27 09:44:35 kmura Exp $
# ===========================================================
#   Makefile for building Python module with Boost.Python
#      
# ===========================================================
include ../../../config/config.gmk

# python module name
MODULE := demo_wp#

# location of exposured library
G4PY_EXLIB_ROOT    := ./g4app
G4PY_EXLIB_INCDIR  := $(G4PY_EXLIB_ROOT)/include
G4PY_EXLIB_LIBDIR  := $(G4PY_EXLIB_ROOT)/tmp/$(G4SYSTEM)/demo_wp
G4PY_EXLIB_LIBNAME := demo_wp

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/module.gmk

.PHONY: all mylib cleanall

all : cleanall mylib module

mylib :
	@cd $(G4PY_EXLIB_ROOT) && $(MAKE)

cleanall: clean
	@cd $(G4PY_EXLIB_ROOT) && $(MAKE) clean
