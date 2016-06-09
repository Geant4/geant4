# $Id: GNUmakefile,v 1.1 2007/11/16 14:28:34 kmura Exp $
# $Name: geant4-09-01 $
# ===========================================================
#   Makefile for building G4 MPI aplication
# ===========================================================
name := exMPI02#
G4TARGET := $(name)
G4EXLIB := true

G4WORKDIR := .

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
include ../mpi_interface/G4MPI.gmk

LDLIBS1 += -lG4UImpi

ifdef ROOTSYS
  CPPFLAGS += -I$(shell $(ROOTSYS)/bin/root-config --incdir)
  LDLIBS   += $(shell $(ROOTSYS)/bin/root-config --libs)
endif

ifeq ($(G4SYSTEM),Linux-g++)
  LDLIBS   += -Wl,-rpath,$(shell $(ROOTSYS)/bin/root-config --libdir)
endif

