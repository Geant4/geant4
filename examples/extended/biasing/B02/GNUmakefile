# $Id: GNUmakefile,v 1.7 2002-11-19 15:03:58 mdressel Exp $
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


# set this variable to your Anaphe instalation.
ANAPHE_SCRIPTS_PATH:=share/LHCXX/5.0.1/scripts

CPPFLAGS +=
#LDFLAGS += -L/usr/local/lib

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += `$(ANAPHETOP)/$(ANAPHE_SCRIPTS_PATH)/aida-config -i`
LDFLAGS += `$(ANAPHETOP)/$(ANAPHE_SCRIPTS_PATH)/aida-config -l`
