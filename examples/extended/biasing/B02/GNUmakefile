# $Id: GNUmakefile,v 1.6 2002-11-18 13:23:01 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------


# this GNUmakefile shows how to use the example together with Anaphe
# you may also chose another AIDA compliant analysis package.

name := exampleB02
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin


# set this variable to your Anaphe instalation.
ANAPHEPATH:=/afs/cern.ch/sw/lhcxx/share/LHCXX/5.0.1/scripts

CPPFLAGS +=

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += `$(ANAPHEPATH)/aida-config -i`
LDFLAGS += `$(ANAPHEPATH)/aida-config -l`
