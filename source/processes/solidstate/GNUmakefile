# $Id: GNUmakefile $
# ------------------------------------------------------------------
# GNUmakefile for solidstate library.  Enrico Bagli, 6/4/2017.
# ------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4solidstate

SUBDIRS  = phonon
SUBDIRS += channeling

SUBLIBS  = G4phonon
SUBLIBS += G4channeling

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
