# $Id: GNUmakefile,v 1.2 1999-05-31 17:38:38 stesting Exp $
# ------------------------------------------------------------------
# GNUmakefile for electromagnetic library.  Gabriele Cosmo, 18/9/96.
# ------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4electromagnetic

SUBDIRS = muons standard utils xrays lowenergy
SUBLIBS = G4muons G4emstandard G4emutils G4xrays G4emlowenergy

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
