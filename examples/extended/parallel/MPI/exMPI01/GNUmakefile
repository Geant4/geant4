# $Id: GNUmakefile,v 1.2 2010/01/11 16:53:51 gcosmo Exp $
# $Name: geant4-09-04-beta-01 $
# ===========================================================
#   Makefile for building G4 MPI aplication
# ===========================================================

name := exMPI01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
include ../mpi_interface/G4MPI.gmk

LDLIBS1 += -lG4UImpi

