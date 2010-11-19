# $Id: GNUmakefile,v 1.9 2010-11-19 08:44:20 gcosmo Exp $
# ------------------------------------------------------------------
# GNUmakefile for electromagnetic library.  Gabriele Cosmo, 18/9/96.
# ------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4electromagnetic

SUBDIRS = muons standard utils xrays lowenergy highenergy adjoint polarisation pii
SUBLIBS = G4muons G4emstandard G4emutils G4xrays G4emlowenergy G4emhighenergy G4emadjoint G4empolar G4empii

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
