# $Id: GNUmakefile,v 1.1 2007/11/16 14:16:46 kmura Exp $
# $Name: geant4-09-02 $
# ===========================================================
#   Makefile for building G4 MPI aplication
# ===========================================================
name := exMPI01#
G4TARGET := $(name)
G4EXLIB := true

G4WORKDIR := .

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
include ../mpi_interface/G4MPI.gmk

LDLIBS1 += -lG4UImpi

