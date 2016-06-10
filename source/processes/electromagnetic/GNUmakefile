# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
# ------------------------------------------------------------------
# GNUmakefile for electromagnetic library.  Gabriele Cosmo, 18/9/96.
# ------------------------------------------------------------------

MAKEFLAGS= --no-print-directory

name := G4electromagnetic

SUBDIRS  = muons standard utils xrays lowenergy
SUBDIRS += pii highenergy adjoint polarisation
SUBDIRS += dna/processes dna/models dna/utils dna/management
SUBDIRS += dna/molecules/management dna/molecules/types

SUBLIBS  = G4muons G4emstandard G4emutils G4xrays G4emlowenergy
SUBLIBS += G4empii G4emhighenergy G4emadjoint G4empolar
SUBLIBS += G4emdna-processes G4emdna-models G4emdna-utils G4emdna-man
SUBLIBS += G4emdna-molman G4emdna-moltypes

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/globlib.gmk
