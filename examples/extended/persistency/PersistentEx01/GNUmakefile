# $Id: GNUmakefile,v 1.5 2001/03/15 22:56:27 morita Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------
# You must have unique G4EXAMPLE_FDID assigned by your local 
# manager of Objectivity/DB lock server.
#

name := PersistentEx01
G4ODBMS := true
G4TARGET := $(name)
G4EXLIB := true
PWD := $(shell (dirname `pwd`/dummy | sed -e 's/ //g'))

ifndef G4INSTALL
  G4INSTALL = $(PWD)/../../../..
endif

.PHONY: all

all: lib bin

include $(G4INSTALL)/config/architecture.gmk

G4OODDLX_BOOT := $(G4EXAMPLE_BOOT)

  CPPFLAGS += \
              -I$(PWD)/include \
              -I$(G4TMPDIR)

include $(G4INSTALL)/config/binmake.gmk

