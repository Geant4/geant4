# $Id: GNUmakefile,v 1.1 2006-02-27 10:05:24 kmura Exp $
# ===========================================================
#   Makefile for building Python module with Boost.Python
#      
# ===========================================================
include ../../config/config.gmk

# python module name
MODULE := gtest01#

# location of exposured library
G4PY_EXLIB_ROOT    := ./g4app
G4PY_EXLIB_INCDIR  := $(G4PY_EXLIB_ROOT)/include
G4PY_EXLIB_LIBDIR  := $(G4PY_EXLIB_ROOT)/tmp/$(G4SYSTEM)/Q00
G4PY_EXLIB_LIBNAME := Q00

include $(G4PY_INSTALL)/config/g4py.gmk
include $(G4PY_INSTALL)/config/module.gmk

.PHONY: all mylib cleanall

all : cleanall mylib module

mylib :
	@cd $(G4PY_EXLIB_ROOT) && $(MAKE)

cleanall: clean
	@cd $(G4PY_EXLIB_ROOT) && $(MAKE) clean
