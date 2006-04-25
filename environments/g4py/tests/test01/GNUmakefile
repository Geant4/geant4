# $Id: GNUmakefile,v 1.2 2006-04-25 07:54:28 kmura Exp $
# ===========================================================
#   Makefile for building Python module with Boost.Python
#      
# ===========================================================
include ../../config/config.gmk

# python module name
MODULE := test01#

# location of exposured library
G4PY_EXLIB_ROOT    := ./mylib
G4PY_EXLIB_INCDIR  := $(G4PY_EXLIB_ROOT)
G4PY_EXLIB_LIBDIR  := $(G4PY_EXLIB_ROOT)
G4PY_EXLIB_LIBNAME := XXX

COMPILE_WITH_BOOST := true
include $(G4PY_INSTALL)/config/module.gmk

.PHONY: all mylib cleanall

all : cleanall mylib module

mylib :
	@cd $(G4PY_EXLIB_ROOT) && $(MAKE)

cleanall: clean
	@cd $(G4PY_EXLIB_ROOT) && $(MAKE) clean

