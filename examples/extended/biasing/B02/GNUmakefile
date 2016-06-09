# $Id: GNUmakefile,v 1.9 2004/06/18 15:59:39 daquinog Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------


# this GNUmakefile shows how to use the example together with Anaphe
# you may also chose another AIDA compliant analysis package.
#
# The environment variable ANAPHETOP has to be setup. See ../README
#

name := exampleB02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
 CPPFLAGS += `aida-config --include`
 LDFLAGS += `aida-config --lib`
endif 
