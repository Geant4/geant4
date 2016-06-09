# $Id: GNUmakefile,v 1.8 2003/06/20 15:32:09 dressel Exp $
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



CPPFLAGS +=
#LDFLAGS += -L/usr/local/lib

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += `aida-config -i`
LDFLAGS += `aida-config -l`
