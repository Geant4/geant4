# $Id: GNUmakefile,v 1.1 1999-01-07 16:11:01 gunter Exp $
# ------------------------------------------------------------------
# GNUmakefile for electromagnetic library.  Gabriele Cosmo, 18/9/96.
# ------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4electromagnetic

SUBDIRS = muons standard utils xrays
SUBLIBS = G4muons G4emstandard G4emutils G4xrays

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
