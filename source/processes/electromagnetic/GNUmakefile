# $Id: GNUmakefile,v 1.3 2001-03-05 10:52:05 maire Exp $
# ------------------------------------------------------------------
# GNUmakefile for electromagnetic library.  Gabriele Cosmo, 18/9/96.
# ------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4electromagnetic

SUBDIRS = muons standard utils xrays lowenergy integral
SUBLIBS = G4muons G4emstandard G4emutils G4xrays G4emlowenergy G4integral

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
