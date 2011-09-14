# $Id: GNUmakefile,v 1.9 2010-11-19 08:44:20 gcosmo Exp $
# ------------------------------------------------------------------
# GNUmakefile for electromagnetic library.  Gabriele Cosmo, 18/9/96.
# ------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4electromagnetic

SUBDIRS  = muons standard utils xrays lowenergy
SUBDIRS += pii highenergy adjoint polarisation
SUBDIRS += dna/processes dna/models dna/utils
SUBLIBS  = G4muons G4emstandard G4emutils G4xrays G4emlowenergy
SUBLIBS += G4empii G4emhighenergy G4emadjoint G4empolar
SUBLIBS += G4emdna-processes G4emdna-models G4emdna-utils

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
